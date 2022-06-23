#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(max);

# -----------------------------------------------------------------------------
# VARIABLES

use vars qw ($USAGE $IN $OUT_P $OUT_M $OUT_S $OUT_E $OUT_X $bam $temp $mate $mode $silent);
use vars qw ($out_paired_freq $out_mate_freq $out_single_freq $out_edist $out_split);
use vars qw (%mapping_freq_paired %mapping_freq_mate %mapping_freq_single);
use vars qw (%edist %XL_freq %XA_SU_freq);
use vars qw ($NM $NH $XF $XB $XQ $XL $XA $XR);
use vars qw ($qname $flag $rname $pos $mapq $cigar $mname $mpos $isize $seq $qual @tags);

my ($unique_paired, $multiple_paired, $unique_mate, $multiple_mate, $unique_single, $multiple_single, $not_split, $split, $nr_mappings, $alter_XA) = (0)x10;
my $NM_out = "WARNING: No NM tag -> no edistance distribution will be calculated!!!\n";
my $XA_out = "WARNING: There are mappings with one end mapped as single mate and both ends are mapped as pair (at least one end is split read). Amount of unique mapped reads might be overestimated.\n";

# -----------------------------------------------------------------------------
# OPTIONS

$USAGE = << "USE";

  usage:  perl count_mapped_reads.pl --bam <file> [-b] [--silent] [-p <file>] [-m <file>] [-s <file>] [-e <file>] [-x <file>]

  options:  -bam    Path to input file (default: stdin)
            -b      if bisulfite mapping
            -p      Path to output file for multiple hits frequencies paired end (default: stdout)
            -m      Path to output file for multiple hits frequencies one mate (default: stdout)
            -s      Path to output file for multiple hits frequencies single ends (default: stdout)
            -e      Path to output file for e-distance distribution (default: stdout)
            -x      Path to output file for split fragment distribution (default: stdout)

  input:    (1a) bam file
                perl count_mapped_reads.pl -bam xxx.bam
                
                or
                
            (1b) sam format input by stdin generated with samtools view
                samtools view xxx.bam | perl count_mapped_reads.pl

  output:   stdout
            (0) total amount of mappings
            (1) mapped reads:
                number of mapped reads, seperated by both ends of a pair
                are mapped once (unique paired ends) or multiple times
                (multiple paired ends) and only one end of a pair mapped
                once (unique one mate) or multiple times (multiple
                one mate) and a single end read mapped once (unique
                single end) or multiple times (multiple single end)
            (2) sum of all mapped reads and if paired ends sum of all
                mapped fragments
            (3) if split reads are present, amount of split and not split
                mates and amount of NH flags set to 1
            
            stdout or provided files
            (4) frequencies of the multiple hits seperated by paired end
                reads and mate pairs and single end reads
                Format: nr of hits<tab>frequency
            (5) frequencies of the e-distances of the reads
                e-distance of ALL (not split) reads is count!
                Format: e-distance<tab>frequency
            (6) frequency of number of split fragments
                Format: nr split fragments<tab>frequency
                
            if option silent
                is set, (0)-(3) are given as tab-seperated numbers and
                (4)-(6) are print to stdout or provied files:
                total amount of mappings<tab>mapped unique paired ends<tab>
                mapped multiple paired ends<tab>mapped unique one mate<tab>
                mapped multiple one mate<tab>mapped unique single end<tab>
                mapped multiple single end<tab>sum of mapped reads<tab>
                sum of mapped mates<tab>amount of split mates<tab>
                amount of not split mates<tab>amount of NH flags set to 1
                
USE

unless (GetOptions(
    "bam=s"   => \$bam,
    "b"       => \$mode,
    "p=s"     => \$out_paired_freq,
    "m=s"     => \$out_mate_freq,
    "s=s"     => \$out_single_freq,
    "e=s"     => \$out_edist,
    "x=s"     => \$out_split,
    "silent"  => \$silent
)){
    printf STDERR $USAGE;
    exit -1;
}

# -----------------------------------------------------------------------------
# MAIN

## if input okay, open stream
if (!defined $bam && -t STDIN){
    printf STDERR $USAGE;
    exit -1;
}
elsif (defined $bam){
    $temp = $bam;
    $temp =~ s/.bam//;
    print "-----------$bam-----------\n\n";
    open($IN, "samtools view $bam | ") or die "Error in input\n";
}
else{
    $IN = \*STDIN;
}

while(<$IN>){
    chomp;
    
    ## get sam information
    ($qname, $flag, undef $rname, undef $pos, undef $mapq, $cigar, undef $mname, undef $mpos, undef $isize, undef $seq, undef $qual, @tags) = split(/\t/,$_);
    
    if ($qname =~ m/^@/){
        next;
    }
    
    ($NM, $NH, $XF, $XB, $XQ, $XL, $XA, $XR) = (undef)x8;
    
    ## get flag information
    foreach (@tags){
        my $tag = substr($_, 0, 5, "");
        if ($tag eq "NM:i:"){
            $NM = $_;           #edist
        }
        elsif ($tag eq "NH:i:"){
            $NH = $_;           #freq
        }
        elsif ($tag eq "XF:i:"){
            $XF = $_;           #nr of bisulfite mismatches
        }
        elsif ($tag eq "XB:Z:"){
            $XB = $_;           #bisulfite mapping
        }
        elsif ($tag eq "XQ:i:"){
            $XQ = $_;           #split read
        }
        elsif ($tag eq "XL:i:"){
            $XL = $_;           #nr of split-fragments
        }
        elsif ($tag eq "XA:Z:"){
            $XA = $_;           #
        }
        elsif ($tag eq "XR:i:"){
            $XR = $_;           #file is realigned
        }
    }
    
    ## check for dependencies
    if (defined $XR){
        die "NO realigned files!!!\n";
    }
    
    if (defined $mode){
        if (!defined $XB){
            die "No XB flag, although -b option is set: $_\n";
        }
        if (defined $XQ){
            die "XQ flag, although -b option is set: $_\n";
        }
        if (!defined $NM){
            die "No XB flag, although -b option is set: $_\n";
        }
    }
    else{
        if (defined $XB){
            die "XB flag, although -b option is not set: $_\n";
        }
        if (!defined $XL && defined $XQ){
            die "Split read without XL flag: $_\n";
        }
    }
    
    if (!defined $NH){
        die "No NH tag: $_\n";
    }
    
    if (!defined $NM){
        if ($NM_out ne 'done'){
            print STDERR $NM_out;
            $NM_out = "done";
        }
    }
    
    if (!defined $XA){
        die "No XA tag: $_\n";
    }
    
    if ($XA =~ m/^[SU]$/ && $NH != 1){
        if ($XA_out ne 'done'){
            print STDERR $XA_out;
            $XA_out = "done";
        }
        if ($NH != 1){
            $alter_XA++;
        }
    }

    ## set mate information
    if ($flag & 0x4){           #read unmapped
        die "Unmapped read: $_\n";
    }
    elsif ($flag & 0x1){
        if ($flag & 0x8){
            $mate = 1;          #only one mate mapped, bislufite: possibly both mates are mapped independently as single mates
        }
        elsif (($flag & 0x40) && ($flag & 0x80)){
            die "Read is first and second in pair: $_\n";
        }
        elsif ($flag & 0x40){
            $mate = 2;          #both ends mapped, first in pair
        }
        elsif ($flag & 0x80){
            $mate = 2;          #both ends mapped, second in pair
        }
        else{
            die "paired-end and not 0x40 or 0x80: $_\n";
        }
    }
    elsif (($flag & 0x40) || ($flag & 0x80)){
        die "not paired-end but first in pair or second in pair is set: $_\n";
    }
    else{
        $mate = 0;              #single ends
    }
    
    ## count mapping frequency
    if (defined $mode || !defined $XQ || $XQ == 0){
        if ($mate == 2){            
            if ($XA =~ m/^[SU]$/){
                if (exists $XA_SU_freq{$NH}){
                    $XA_SU_freq{$NH}++;
                }
                else{
                    $XA_SU_freq{$NH} = 1;
                }
            }
            else {            
                if (exists $mapping_freq_paired{$NH}){  #paired end
                    $mapping_freq_paired{$NH}++;
                }
                else{
                    $mapping_freq_paired{$NH} = 1;
                }
            }
        }
        elsif ($mate == 1){                         #mate
            if (exists $mapping_freq_mate{$NH}){
                $mapping_freq_mate{$NH}++;
            }
            else{
                $mapping_freq_mate{$NH} = 1;
            }
        }
        elsif ($mate == 0){                         #single ends
            if (exists $mapping_freq_single{$NH}){
                $mapping_freq_single{$NH}++;
            }
            else{
                $mapping_freq_single{$NH} = 1;
            }
        }
        
        ## count total amount of mappings
        $nr_mappings++;
    }
    
    ## count edists of ALL (not split) reads
    if (defined $mode || !defined $XQ){
        if (exists $edist{$NM}){
            $edist{$NM}++;
        }
        else{
            $edist{$NM} = 1;
        }
    }
    
    ## count number of split-fragment frequency
    if (defined $XQ && $XQ == 0){
        if (exists $XL_freq{$XL}){
            $XL_freq{$XL}++;
        }
        else{
            $XL_freq{$XL} = 1;
        }
    }
}

## close stream
if (defined $bam){
    close($IN);
}

foreach (keys %XA_SU_freq){
    $mapping_freq_paired{1} = $mapping_freq_paired{1}+$XA_SU_freq{$_}/$_;
}

## count number of uniquely mapped paired ends
if (exists $mapping_freq_paired{1}){
    $unique_paired = ($mapping_freq_paired{1}/2);
    delete ($mapping_freq_paired{1});
}
else{
    $unique_paired = 0;
}
## count number of multiple mapped paired ends
foreach (keys %mapping_freq_paired){
    $multiple_paired = ($multiple_paired+($mapping_freq_paired{$_}/(2*$_)));
}

## count number of uniquely mapped mates
if (exists $mapping_freq_mate{1}){
    $unique_mate = $mapping_freq_mate{1};
    delete ($mapping_freq_mate{1});
}
else{
    $unique_mate = 0;
}
## count number of multiple mapped mates
foreach (keys %mapping_freq_mate){
    $multiple_mate = ($multiple_mate+($mapping_freq_mate{$_}/$_));
}

## count number of uniquely mapped single
if (exists $mapping_freq_single{1}){
    $unique_single = $mapping_freq_single{1};
    delete ($mapping_freq_single{1});
}
else{
    $unique_single = 0;
}
## count number of multiple mapped single
foreach (keys %mapping_freq_single){
    $multiple_single = ($multiple_single+($mapping_freq_single{$_}/$_));
}

## calculate sums of mapped reads/fragments
my $sum_reads = $unique_paired+$multiple_paired+$unique_mate+$multiple_mate+$unique_single+$multiple_single;
my $sum_frag  = (2*$unique_paired)+(2*$multiple_paired)+$unique_mate+$multiple_mate+$unique_single+$multiple_single;

## print 'important' numbers, if silent
if (defined $silent){
    print "$nr_mappings\t$unique_paired\t$multiple_paired\t$unique_mate\t$multiple_mate\t$unique_single\t$multiple_single\t$sum_reads\t$sum_frag\t$split\t$not_split\t$alter_XA\n";
}

## print output: number of mapped reads
if (!defined $silent){
    print "-----------------------------\n--------   SUMMARY   --------\n-----------------------------\n";
    print "\ntotal amount of MAPPINGS:\t$nr_mappings\n";
    print "\nMAPPED READS:\nunique paired end:\t$unique_paired\nmultiple paired end:\t$multiple_paired\nunique one mate:\t$unique_mate\nmultiple one mate:\t$multiple_mate\nunique single end:\t$unique_single\nmultiple single end:\t$multiple_single\n\n";
    print "sum of mapped reads:\t$sum_reads\n";
    if ($sum_reads ne $sum_frag){
        print "sum of mapped mates:\t$sum_frag\n";
    }
    if (!defined $mode){
        print "\namount of split mates(reads):\t\t$split\n";
        print "amount of not split mates(reads):\t$not_split\n";
    }
    if ($alter_XA != 0){
        print "\namount of NH flags set to 1:\t$alter_XA\n";
    }
    print "\n-----------------------------\n\n";
}

## print output: mapping frequency for paired end reads
if (defined $out_paired_freq){
    open($OUT_P, ">$out_paired_freq") or die "Error in output paired end frequencies\n";
}
else{
    $OUT_P = \*STDOUT;
    if (%mapping_freq_paired || $unique_paired){
        print $OUT_P "paired_frequencies:\nhits\tfrequency\n1\t$unique_paired\n";
    }
}
foreach (sort{$a <=> $b} keys %mapping_freq_paired){
    my $nr = ($mapping_freq_paired{$_}/(2*$_));
    print $OUT_P "$_\t$nr\n";
}
if (defined $out_paired_freq){
    close($OUT_P);
}

## print output: mapping frequency for mates
if (defined $out_mate_freq){
    open($OUT_M, ">$out_mate_freq") or die "Error in output mate frequencies\n";
}
else{
    $OUT_M = \*STDOUT;
    if (%mapping_freq_mate || $unique_mate){
        print $OUT_M "\nmate_frequencies:\nhits\tfrequency\n1\t$unique_mate\n";
    }
}
foreach (sort{$a <=> $b} keys %mapping_freq_mate){
    my $nr = ($mapping_freq_mate{$_}/$_);
    print $OUT_M "$_\t$nr\n";
}
if (defined $out_mate_freq){
    close($OUT_M);
}

## print output: mapping frequency for single ends
if (defined $out_single_freq){
    open($OUT_S, ">$out_single_freq") or die "Error in output single frequencies\n";
}
else{
    $OUT_S = \*STDOUT;
    if (%mapping_freq_single || $unique_single){
        print $OUT_S "\nsingle_frequencies:\nhits\tfrequency\n1\t$unique_single\n";
    }
}
foreach (sort{$a <=> $b} keys %mapping_freq_single){
    my $nr = ($mapping_freq_single{$_}/$_);
    print $OUT_S "$_\t$nr\n";
}
if (defined $out_single_freq){
    close($OUT_S);
}

## print output: e-distance for ALL (not split) reads
if (defined $out_edist){
    open($OUT_E, ">$out_edist") or die "Error in output e-distance\n";
}
else{
    $OUT_E = \*STDOUT;
    if (%edist){
        print $OUT_E "\ne_distance of ALL (not split) mappings:\nedist\tfrequency\n";
    }
}
foreach (sort{$a <=> $b} keys %edist){
    print $OUT_E "$_\t$edist{$_}\n";
}
if (defined $out_edist){
    close($OUT_E);
}

## print output: number of split fragments
if (defined $out_split){
    open($OUT_X, ">$out_split") or die "Error in output nr_fragments\n";
}
else{
    $OUT_X = \*STDOUT;
    if (%XL_freq){
        print $OUT_X "\nnumber of split reads fragments:\nnr_fragments\tfrequency\n";
    }
}
foreach (sort{$a <=> $b} keys %XL_freq){
    print $OUT_X "$_\t$XL_freq{$_}\n";
    $split     = $split+$XL_freq{$_};
}
if (defined $out_split){
    close($OUT_X);
}
$not_split = $sum_frag-$split;

