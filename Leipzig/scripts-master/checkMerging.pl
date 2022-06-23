#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max min);
use vars qw/$USAGE %FILES %KMER $K %PREV %ENTRY $CLIPTHR %STATS $ADAPTER $ADAPTER2 $MINOVERLAP $INDEL/;

$ADAPTER="CTGTCTCTTATACACATCTCCGAGCCCACGAGACTAAGGCGAATCTCGACAGCCGTCTTCTGCTTG";
$ADAPTER2="CTGTCTCTTATACACATCTGACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT";
$CLIPTHR = 10;
$INDEL = -2;
$MINOVERLAP = 30;


$USAGE = << "USE";
  usage:  perl processMerging.pl -i <file> -j <file> -k <file> -l <file> -m <file> -x <file> -o <prefix>

  options:  -i     unclipped read 1 fastQ file
            -j     unclipped read 2 fastQ file
            -k     clipped read 1 fastQ file
            -l     clipped read 2 fastQ file
            -m     clipped and merged reads fastQ file
            -x     k-mer count file
            -o     prefix of output files

  input:    see README

  output:   see README
 
  version:  0.1

 written by Christian Otto, Bioinformatics Leipzig, Dec 2012
USE

unless (GetOptions(
	    "i=s"  => \$FILES{raw0_f},
            "j=s"  => \$FILES{raw1_f},
	    "k=s"  => \$FILES{clip0_f},
	    "l=s"  => \$FILES{clip1_f},
            "m=s"  => \$FILES{merg_f},
	    "o=s"  => \$FILES{out_f},
	    "x=s"  => \$FILES{kmer_f})){
    printf STDERR $USAGE;
    exit -1;
}

if (!defined $FILES{raw0_f} || !defined $FILES{raw1_f} || !defined $FILES{clip0_f} || 
    !defined $FILES{clip1_f} || !defined $FILES{merg_f} || !defined $FILES{out_f}){
    printf STDERR $USAGE;
    exit -1;
}

sub main {
    ## open input file handles
    $FILES{raw0}  = openFH($FILES{raw0_f}, 1);
    $FILES{raw1}  = openFH($FILES{raw1_f}, 1);
    $FILES{clip0} = openFH($FILES{clip0_f}, 1);
    $FILES{clip1} = openFH($FILES{clip1_f}, 1);
    $FILES{merg}  = openFH($FILES{merg_f}, 1);
    ## open output file handes
    openOutFH($FILES{out_f});
    
    ## read k-mer statistics
    if (defined $FILES{kmer_f}){
	readKmer();
    }
    
    ## init statistics
    initStat(\%STATS);

    ## process input
    while(readNext()){
	if (defined $ENTRY{clip0} && defined $ENTRY{clip1}){
	    processEntry(\%ENTRY, \%STATS);
	}
    }
    ## print statistics
    printStat(\%STATS);

    ## close file handles
    foreach(keys(%FILES)){
	close($_) if !/_f$/;
    }
}

sub readKmer {
    $FILES{kmer}  = openFH($FILES{kmer_f}, 0);
    while(readline($FILES{kmer})){
	chomp;
	my @F = split;
	$KMER{$F[0]} = $F[1];
	if (!defined $K){
	    $K = length($F[0]);
	}
	die "Error in kmer: $_ $K ".length($F[0])."\n" if $K != length($F[0]);
    }
}

sub readNext {

    $PREV{raw0}= readline($FILES{raw0});
    $PREV{raw1}= readline($FILES{raw1});

    $PREV{clip0} = readline($FILES{clip0}) if !defined $PREV{clip0};
    $PREV{clip1} = readline($FILES{clip1}) if !defined $PREV{clip1};
    $PREV{merg}  = readline($FILES{merg})  if !defined $PREV{merg};

    if (!defined $PREV{raw0} || !defined $PREV{raw1} ||
	!defined $PREV{clip0} || !defined $PREV{clip1}){
	die "Error in merged fastQ file: Entry '$ENTRY{merg}->[0]' not yet processed. File not sorted properly?" if defined $PREV{merg};	
	return 0;
    }

    chomp($PREV{raw0}, $PREV{raw1}, $PREV{clip0}, $PREV{clip1});
    chomp($PREV{merg}) if defined $PREV{merg};

    $ENTRY{raw0}  = [split("\t", $PREV{raw0})];
    $ENTRY{raw1}  = [split("\t", $PREV{raw1})];
    
    $ENTRY{clip0} = [split("\t", $PREV{clip0})];
    $ENTRY{clip1} = [split("\t", $PREV{clip1})];
    $ENTRY{merg}  = [split("\t", $PREV{merg})] if defined $PREV{merg};

    #print "Processing:".$ENTRY{raw0}->[0].", next merged:".$ENTRY{merg}->[0]."\n";
    
    ## check description of raw with clipped
    ## whether entry was discarded during clipping
    if ($ENTRY{raw0}->[0] ne $ENTRY{clip0}->[0]){
    	die "Error in parsing\n" if $ENTRY{raw1}->[0] eq $ENTRY{clip1}->[0] || 
    	    $ENTRY{raw0}->[0] eq $ENTRY{merg}->[0];
    	#print $ENTRY{raw0}->[0]." discarded during clipping.\n";

	undef $ENTRY{clip0};
	undef $ENTRY{clip1};
	undef $ENTRY{merg};
    	return 1;
    }
    die "Error in parsing\n" if $ENTRY{raw1}->[0] ne $ENTRY{clip1}->[0];
        
    ## next entry is merged
    if (defined $PREV{merg} && 
	$ENTRY{raw0}->[0] eq $ENTRY{merg}->[0]){
    	undef $PREV{merg};
    }
    else {
	undef $ENTRY{merg};
    }
    undef $PREV{clip0};
    undef $PREV{clip1};

    return 1;
}

sub openFH {
    my ($f, $fq) = @_;
    my $fh;
    if ($f =~ /\.gz$/){
	if ($fq){
	    open($fh, "gunzip -c $f | paste - - - - |") or die;
	}
	else {
	    open($fh, "gunzip -c $f |") or die;
	}
    }
    else {
	if ($fq){
	    open($fh, "cat $f | paste - - - - |") or die;
	}
	else {
	    open($fh, "< $f") or die;
	}
    }
    return $fh;
}

sub openOutFH {
    my ($f) = @_;

    open($FILES{out0}, "| gzip -c > ".$FILES{out_f}.".R1.fastq.gz") or die;
    open($FILES{out1}, "| gzip -c > ".$FILES{out_f}.".R2.fastq.gz") or die;
    open($FILES{outm}, "| gzip -c > ".$FILES{out_f}.".M.fastq.gz") or die;
}

sub checkKmer {
    my ($str) = @_;

    return undef if !defined $K || length($str) < $K;

    for (my $i = 0; $i < length($str)-$K+1; $i++){
	my $substr = substr($str, $i, $K);
	if (defined $KMER{$substr}){
	    return $substr."\t".$substr."\t".$KMER{$substr}."\n";
	}
	if (defined $KMER{revcomp($substr)}){
	    return $substr."\t".revcomp($substr)."\t".$KMER{revcomp($substr)}."\n";
	}
    }
    return undef;
}

sub processEntry {
    my ($entry_ref, $stat_ref) = @_;
    my ($raw0, $raw0len, $raw1, $raw1len, $clip0, $clip0len,
	$clip1, $clip1len, $merg, $merglen);

    $raw0 = $entry_ref->{raw0}->[1];
    $raw0len = length($raw0);
    $raw1 = $entry_ref->{raw1}->[1];
    $raw1len = length($raw1);
    $clip0 = $entry_ref->{clip0}->[1];
    $clip0len = length($clip0);
    $clip1 = $entry_ref->{clip1}->[1];
    $clip1len = length($clip1);
    if (defined $entry_ref->{merg}){
	$merg = $entry_ref->{merg}->[1];
	$merglen = length($merg);
    }
    

    $stat_ref->{"all"}++;

    ## case 1: merged
    if (defined $entry_ref->{merg}){
	$stat_ref->{"1"}++;

	## case 1.1: both not clipped
	if ($raw0len == $clip0len && $raw1len == $clip1len){
	    $stat_ref->{"1.1"}++;

	    my ($overlap, $ret);

	    $overlap = substr($raw0, $merglen-$raw1len);
	    $ret = checkKmer($overlap);

	    ## K-mer file not defined or overlap is not repetitive
	    if (!defined $ret){
		$stat_ref->{"1.1.1"}++;
		
		## report merged
		reportMerged($entry_ref);
	    }
	    else {
		$stat_ref->{"1.1.2"}++;
		
		## report as unmerged
		reportUnmerged($entry_ref);
	    }
	    return;
	}
	## case 1.2: both clipped
	elsif ($raw0len != $clip0len && $raw1len != $clip1len){
	    $stat_ref->{"1.2"}++;	
	    
	    ## case 1.2.1: no overhang
	    ##
	    if ($merglen - $clip0len == 0 && $merglen - $clip1len == 0){
		$stat_ref->{"1.2.1"}++;	  

		## simply report as merged
		## (no repeat filter necessary)
		reportMerged($entry_ref);

		return;
	    }
	    my (@DP, $scr, $aln0, $aln1, $acc);
	    my ($i, $seq0, $seq1);

	    ## case 1.2.2: incorrectly shifted merging
	    ##
	    ## - calculate global alignment
	    ##
	    @DP = align_matrix($clip0, revcomp($clip1), "global");
	    ($scr, $aln0, $aln1) = backtrace_matrix($clip0, revcomp($clip1), "global", @DP);
	    $acc = 1 - getEdist($aln0, $aln1, "global")/max($clip0len, $clip1len);

	    if ($acc >= 0.90){
		$stat_ref->{"1.2.2"}++;		
		
		## construct merging based on global alignment
		## (no repeat filter necessary)		
		constructMerging($entry_ref, $aln0, $aln1);
		reportMerged($entry_ref);

		#printEntry($entry_ref);
		#printAlign($aln0, $aln1, "global");
		return;
	    }

	    ## case 1.2.3: incorrectly clipped adapter
	    ##
	    ## - reattach clippings by simply taking raw read
	    ## - calculate global alignment of overlap
	    ##
	    if ($raw0len <= $merglen && $raw1len <= $merglen){

		$seq0 = substr($raw0, $merglen - $raw1len);
		$seq1 = revcomp(substr($raw1, $merglen - $raw0len));
	    
		@DP = align_matrix($seq0, $seq1, "global");
		($scr, $aln0, $aln1) = backtrace_matrix($seq0, $seq1, "global", @DP);
		$acc = 1 - getEdist($aln0, $aln1, "global")/max(length($seq0), length($seq1));
	    		     
		if ($acc >= 0.90){	    
		    $stat_ref->{"1.2.3"}++;	

		    my ($overlap, $ret);

		    $overlap = substr($raw0, $merglen-$raw1len);
		    $ret = checkKmer($overlap);

		    #printEntry($entry_ref);
		    #printAlign($aln0, $aln1, "global");

		    ## K-mer file not defined or overlap is not repetitive
		    if (!defined $ret){
			$stat_ref->{"1.2.3.1"}++;
			
			## report merged
			reportMerged($entry_ref);

			#print "not repetitive:\n";
			#printEntry($entry_ref);
		    }
		    else {
			$stat_ref->{"1.2.3.2"}++;
			
			## report as unmerged
			reportUnmerged($entry_ref);

			#print "repetitive:\n";
			#printEntry($entry_ref);
		    }
		    return;
		}
	    }
	    
	    ## case 1.2.4:
	    ## - report as unmerged
	    ##
	    $stat_ref->{"1.2.4"}++;
	    
	    reportUnmerged($entry_ref);
	    return;
	}
	## case 1.3: either of them clipped
	else {
	    $stat_ref->{"1.3"}++;	  
  
	    ## case 1.3.1: incorrectly clipped adapter
	    ##
	    ## - reattach clippings by simply taking raw read
	    ## - recalculate overlap alignment
	    ##
	    my (@DP, $scr, $aln0, $aln1, $acc);	    
	    my ($i, $len, $seq0, $seq1);
	    
	    if ($raw0len <= $merglen && $raw1len <= $merglen){
	    
	    	$seq0 = substr($raw0, $merglen - $raw1len);
	    	$seq1 = revcomp(substr($raw1, $merglen - $raw0len));
	    
		@DP = align_matrix($seq0, $seq1, "global");
		($scr, $aln0, $aln1) = backtrace_matrix($seq0, $seq1, "global", @DP);
		$acc = 1 - getEdist($aln0, $aln1, "global")/max(length($seq0), length($seq1));
	    
		if ($acc >= 0.90){
		    $stat_ref->{"1.3.1"}++;

		    my ($overlap, $ret);

		    $overlap = substr($raw0, $merglen-$raw1len);
		    $ret = checkKmer($overlap);

		    #printEntry($entry_ref);
		    #printAlign($aln0, $aln1, "global");

		    ## K-mer file not defined or overlap is not repetitive
		    if (!defined $ret){
			$stat_ref->{"1.3.1.1"}++;
			
			## report merged
			reportMerged($entry_ref);

			#print "not repetitive:\n";
			#printEntry($entry_ref);
		    }
		    else {
			$stat_ref->{"1.3.1.2"}++;
			
			## report as unmerged
			reportUnmerged($entry_ref);

			#print "repetitive:\n";
			#printEntry($entry_ref);
		    }
		    return;
		}
	    }


	    ## case 1.3.2: otherwise
	    ## - report as unmerged
	    ##
	    $stat_ref->{"1.3.2"}++;

	    ## report as unmerged
	    reportUnmerged($entry_ref); 
	    return;
	}
    }
    ## case 2: unmerged
    else {
	$stat_ref->{"2"}++;
	## case 2.1: both not clipped or only short clippings
	##
	## - if any, only short (probably false positive) clippings
	## - report as correctly not merged
	##
	if ($raw0len - $clip0len < $CLIPTHR &&
	    $raw1len - $clip1len < $CLIPTHR){
	    
	    $stat_ref->{"2.1"}++;

	    ## report as unmerged
	    reportUnmerged($entry_ref); 
	    return;
	}
	## case 2.2: both clipped (long clipping)
	elsif ($raw0len - $clip0len >= $CLIPTHR &&
	       $raw1len - $clip1len >= $CLIPTHR){

	    $stat_ref->{"2.2"}++;

	    ## case 2.2.1: clipped sequences too short
	    ##
	    ## - overlap was not reported during merging simply
	    ##   due to minOverlap parameter
	    ## - otherwise, it would have be characterized as case 1.2.1
	    ## - calculate global alignment with very high strictness	 
	    ##
	    if ($clip0len < $MINOVERLAP || $clip1len < $MINOVERLAP){
		
		my (@DP, $scr, $aln0, $aln1, $acc);
		@DP = align_matrix($clip0, revcomp($clip1), "global");
		($scr, $aln0, $aln1) = backtrace_matrix($clip0, revcomp($clip1), "global", @DP);
		$acc = 1 - getEdist($aln0, $aln1, "global")/max($clip0len, $clip1len);
		if ($acc >= 0.95){
		    $stat_ref->{"2.2.1"}++;
		    
		    ## construct merging based on global alignment
		    ## (no repeat filter necessary)
		    constructMerging($entry_ref, $aln0, $aln1);
		    reportMerged($entry_ref);
		    
		    #printEntry($entry_ref);
		    #printAlign($aln0, $aln1, "global");
		    return;
		}
	    }

	    ## case 2.2.2: otherwise
	    ##
	    ## - mostly at least one of the reads is low-quality and hence 
	    ##   merging is not possible with sufficient overlap length
	    ## - alternative: incorrectly clipped adapters (but less likely)
	    ## - report as correctly not merged
	    ##
	    $stat_ref->{"2.2.2"}++;

	    ## report as unmerged
	    reportUnmerged($entry_ref);

	    #printEntry($entry_ref);
	    #printAlign($aln0, $aln1, "global");
	    return;
	}
	## case 2.3: either of them clipped (long clipping)
	else {
	    $stat_ref->{"2.3"}++;

	    my (@DP, $scr, $seq0, $seq1, $aln0, $aln1, $acc);

	    ## case 2.3.1: incorrectly not clipped adapter
	    ##
	    ## - clip not clipped sequence to same length
	    ##   as clipped sequence 
	    ## - calculate global alignment
	    ##
	    if ($raw0len - $clip0len >= $CLIPTHR){
		$seq0 = $clip0;
		$seq1 = revcomp(substr($clip1, 0, $clip0len))
	    }
	    else {
		$seq0 = substr($clip0, 0, $clip1len);
		$seq1 = revcomp($clip1);

	    }
	    @DP = align_matrix($seq0, $seq1, "global");
	    ($scr, $aln0, $aln1) = backtrace_matrix($seq0, $seq1, "global", @DP);
	    $acc = 1 - getEdist($aln0, $aln1, "global")/length($seq0);

	    if ($acc >= 0.95) {
		$stat_ref->{"2.3.1"}++;
		
		## construct merging based on global alignment
		## (no repeat filter necessary)	
		constructMerging($entry_ref, $aln0, $aln1);
		reportMerged($entry_ref);
		
		#printEntry($entry_ref);
		#printAlign($aln0, $aln1, "global");
		return;
	    }	    

	    ## case 2.3.2: otherwise
	    ##
	    ## - mostly at least one of the reads is low-quality and hence 
	    ##   merging is not possible with sufficient overlap length
	    ## - alternative: incorrectly clipped adapters (but less likely)
	    ## - report as correctly not merged
	    ##
	    $stat_ref->{"2.3.2"}++;

	    ## report as unmerged
	    reportUnmerged($entry_ref);

	    #printEntry($entry_ref);
	    #printAlign($aln0, $aln1, "global");
	    return;
	}
    }
}

sub printStat {
    my ($stat_ref) = @_;

    print STDERR "all\t".$stat_ref->{"all"}."\n"; 
    printf STDERR "1\t%d\t%.2f\n", $stat_ref->{"1"}, (100*$stat_ref->{"1"}/$stat_ref->{"all"});  
    printf STDERR "1.1\t%d\t%.2f\n", $stat_ref->{"1.1"}, (100*$stat_ref->{"1.1"}/$stat_ref->{"all"});
    printf STDERR "1.1.1\t%d\t%.2f\n", $stat_ref->{"1.1.1"}, (100*$stat_ref->{"1.1.1"}/$stat_ref->{"all"});
    printf STDERR "1.1.2\t%d\t%.2f\n", $stat_ref->{"1.1.2"}, (100*$stat_ref->{"1.1.2"}/$stat_ref->{"all"});
    printf STDERR "1.2\t%d\t%.2f\n", $stat_ref->{"1.2"}, (100*$stat_ref->{"1.2"}/$stat_ref->{"all"});
    printf STDERR "1.2.1\t%d\t%.2f\n", $stat_ref->{"1.2.1"}, (100*$stat_ref->{"1.2.1"}/$stat_ref->{"all"});
    printf STDERR "1.2.2\t%d\t%.2f\n", $stat_ref->{"1.2.2"}, (100*$stat_ref->{"1.2.2"}/$stat_ref->{"all"});
    printf STDERR "1.2.3\t%d\t%.2f\n", $stat_ref->{"1.2.3"}, (100*$stat_ref->{"1.2.3"}/$stat_ref->{"all"});
    printf STDERR "1.2.3.1\t%d\t%.2f\n", $stat_ref->{"1.2.3.1"}, (100*$stat_ref->{"1.2.3.1"}/$stat_ref->{"all"});
    printf STDERR "1.2.3.2\t%d\t%.2f\n", $stat_ref->{"1.2.3.2"}, (100*$stat_ref->{"1.2.3.2"}/$stat_ref->{"all"});
    printf STDERR "1.2.4\t%d\t%.2f\n", $stat_ref->{"1.2.4"}, (100*$stat_ref->{"1.2.4"}/$stat_ref->{"all"});
    printf STDERR "1.3\t%d\t%.2f\n", $stat_ref->{"1.3"}, (100*$stat_ref->{"1.3"}/$stat_ref->{"all"});
    printf STDERR "1.3.1\t%d\t%.2f\n", $stat_ref->{"1.3.1"}, (100*$stat_ref->{"1.3.1"}/$stat_ref->{"all"});
    printf STDERR "1.3.1.1\t%d\t%.2f\n", $stat_ref->{"1.3.1.1"}, (100*$stat_ref->{"1.3.1.1"}/$stat_ref->{"all"});
    printf STDERR "1.3.1.2\t%d\t%.2f\n", $stat_ref->{"1.3.1.2"}, (100*$stat_ref->{"1.3.1.2"}/$stat_ref->{"all"});
    printf STDERR "1.3.2\t%d\t%.2f\n", $stat_ref->{"1.3.2"}, (100*$stat_ref->{"1.3.2"}/$stat_ref->{"all"});
    printf STDERR "2\t%d\t%.2f\n", $stat_ref->{"2"}, (100*$stat_ref->{"2"}/$stat_ref->{"all"});  
    printf STDERR "2.1\t%d\t%.2f\n", $stat_ref->{"2.1"}, (100*$stat_ref->{"2.1"}/$stat_ref->{"all"});
    printf STDERR "2.2\t%d\t%.2f\n", $stat_ref->{"2.2"}, (100*$stat_ref->{"2.2"}/$stat_ref->{"all"});
    printf STDERR "2.2.1\t%d\t%.2f\n", $stat_ref->{"2.2.1"}, (100*$stat_ref->{"2.2.1"}/$stat_ref->{"all"});
    printf STDERR "2.2.2\t%d\t%.2f\n", $stat_ref->{"2.2.2"}, (100*$stat_ref->{"2.2.2"}/$stat_ref->{"all"});
    printf STDERR "2.3\t%d\t%.2f\n", $stat_ref->{"2.3"}, (100*$stat_ref->{"2.3"}/$stat_ref->{"all"});
    printf STDERR "2.3.1\t%d\t%.2f\n", $stat_ref->{"2.3.1"}, (100*$stat_ref->{"2.3.1"}/$stat_ref->{"all"});
    printf STDERR "2.3.2\t%d\t%.2f\n", $stat_ref->{"2.3.2"}, (100*$stat_ref->{"2.3.2"}/$stat_ref->{"all"});
}

sub initStat {
    my ($stat_ref) = @_;
    $stat_ref->{"all"} = 0;
    $stat_ref->{"1"} = 0;
    $stat_ref->{"1.1"} = 0;
    $stat_ref->{"1.1.1"} = 0;
    $stat_ref->{"1.1.2"} = 0;
    $stat_ref->{"1.2"} = 0;
    $stat_ref->{"1.2.1"} = 0;
    $stat_ref->{"1.2.2"} = 0;
    $stat_ref->{"1.2.3"} = 0;
    $stat_ref->{"1.2.3.1"} = 0;
    $stat_ref->{"1.2.3.2"} = 0;
    $stat_ref->{"1.2.4"} = 0;
    $stat_ref->{"1.3"} = 0;
    $stat_ref->{"1.3.1"} = 0;
    $stat_ref->{"1.3.1.1"} = 0;
    $stat_ref->{"1.3.1.2"} = 0;
    $stat_ref->{"1.3.2"} = 0;
    $stat_ref->{"2"} = 0;
    $stat_ref->{"2.1"} = 0;
    $stat_ref->{"2.2"} = 0;
    $stat_ref->{"2.2.1"} = 0;
    $stat_ref->{"2.2.2"} = 0;
    $stat_ref->{"2.3"} = 0;
    $stat_ref->{"2.3.1"} = 0;
    $stat_ref->{"2.3.2"} = 0;
}

sub constructMerging {
    my ($entry_ref, $aln0, $aln1) = @_;
    my ($i, @aln0, @aln1, @qual0, @qual1, $seq, $qual);
    die if length($aln0) != length($aln1);

    @aln0 = split(//, $aln0);
    @aln1 = split(//, $aln1);

    $i = 0;
    foreach(@aln0){
	if ($_ ne "-"){
	    push(@qual0, ord(substr($entry_ref->{clip0}->[3], $i++, 1)));
	}
	else {
	    push(@qual0, 0);
	}
    }
    $i = 0;
    foreach(@aln1){
	if ($_ ne "-"){
	    push(@qual1, ord(substr($entry_ref->{clip1}->[3], $i++, 1)));
	}
	else {
	    push(@qual1, 0);
	}
    }
    die if scalar(@aln0) != scalar(@qual0) ||
	scalar(@aln0) != scalar(@qual1);

    #print $entry_ref->{clip0}->[3]."\n";
    #print join(":", @qual0)."\n";
    #print join(":", @aln0)."\n";
    #print join(":", @aln1)."\n";
    #print join(":", @qual1)."\n";
    #print $entry_ref->{clip1}->[3]."\n";
    
    $seq = $qual = "";
    for ($i = 0; $i < length($aln0); $i++){
	if ($aln0[$i] ne "-" && $aln1[$i] ne "-"){
	    if ($qual0[$i] > $qual1[$i]){
		$seq .= $aln0[$i];
		$qual .= chr($qual0[$i]);
	    }
	    else {
		$seq .= $aln1[$i];
		$qual .= chr($qual1[$i]);
	    }
	}
	elsif ($aln0[$i] ne "-"){
	    $seq .= $aln0[$i];
	    $qual .= chr($qual0[$i]);
	}
	elsif ($aln1[$i] ne "-"){
	    $seq .= $aln1[$i];
	    $qual .= chr($qual1[$i]);
	}
	else {
	    die;
	}
    }

    $entry_ref->{merg}->[0] = $entry_ref->{clip0}->[0];
    $entry_ref->{merg}->[1] = $seq;
    $entry_ref->{merg}->[2] = $entry_ref->{clip0}->[2];
    $entry_ref->{merg}->[3] = $qual;
}

sub reportMerged {
    my ($entry_ref) = @_;
    
    print {$FILES{outm}} join("\n", @{$entry_ref->{merg}})."\n";
}

sub reportUnmerged {
    my ($entry_ref) = @_;

    my ($raw0len, $raw1len, $clip0len, $clip1len);

    $raw0len = length($entry_ref->{raw0}->[1]);
    $raw1len = length($entry_ref->{raw1}->[1]);
    $clip0len = length($entry_ref->{clip0}->[1]);
    $clip1len = length($entry_ref->{clip1}->[1]);
   
    ## - report clipped sequences if CLIPTHR is exceeded
    ## - otherwise report raw sequences (to restore FP clippings)    
    if ($raw0len - $clip0len >= $CLIPTHR){
	print {$FILES{out0}} join("\n", @{$entry_ref->{clip0}})."\n";
    }
    else {
	print {$FILES{out0}} join("\n", @{$entry_ref->{raw0}})."\n";
    }
    if ($raw1len - $clip1len >= $CLIPTHR){
	print {$FILES{out1}} join("\n", @{$entry_ref->{clip1}})."\n";
    }
    else {
	print {$FILES{out1}} join("\n", @{$entry_ref->{raw1}})."\n";	
    }
}

sub printEntry {
    my ($entry_ref) = @_;

    my ($raw0, $raw0len, $raw1, $raw1len, $clip0, $clip0len,
	$clip1, $clip1len, $merg, $merglen);

    $raw0 = $entry_ref->{raw0}->[1];
    $raw0len = length($raw0);
    $raw1 = $entry_ref->{raw1}->[1];
    $raw1len = length($raw1);
    $clip0 = $entry_ref->{clip0}->[1];
    $clip0len = length($clip0);
    $clip1 = $entry_ref->{clip1}->[1];
    $clip1len = length($clip1);
    if (defined $entry_ref->{merg}){
	$merg = $entry_ref->{merg}->[1];
	$merglen = length($merg);
    }

    print $entry_ref->{raw0}->[0]."\n";
    if (defined $merg){
	my ($off0, $off1);
	$off0 = max(0, $raw1len - $merglen);
	$off1 = max(0, $merglen - $raw1len);
	print " " x ($off0 + $clip0len), substr($ADAPTER, 0, $raw0len-$clip0len)." ".($raw0len-$clip0len)."\n";
	print " " x $off0, $raw0."\n";
	print " " x $off0, $clip0."\n";
	print " " x $off0, $merg."\n";
	print " " x ($off1+$raw1len-$clip1len), revcomp($clip1)."\n";
	print " " x $off1, revcomp($raw1)."\n";
	print " " x $off1, revcomp(substr($ADAPTER2, 0, $raw1len-$clip1len))." ".($raw1len-$clip1len)."\n\n";
    }
    else {
	print " " x ($clip0len), substr($ADAPTER, 0, $raw0len-$clip0len)." ".($raw0len-$clip0len)."\n";
	print $raw0."\n".$clip0."\n\n";
	print " " x ($raw1len-$clip1len),revcomp($clip1)."\n".revcomp($raw1)."\n";
	print revcomp(substr($ADAPTER2, 0, $raw1len-$clip1len))." ".($raw1len-$clip1len)."\n\n";
    }
}

sub revcomp {
    my ($seq) = @_;
    
    $seq =~ tr/[ACGT]/[TGCA]/;
    $seq = reverse($seq);
    return $seq;
}

sub init_matrix {
    my ($alen, $blen, $type) = @_;
    my ($i, $j, @DP);
    
    if ($type eq "global"){
	for ($i = 0; $i <= $alen; $i++){ $DP[$i][0] = $INDEL * $i }
	for ($j = 0; $j <= $blen; $j++){ $DP[0][$j] = $INDEL * $j }
    }
    elsif ($type eq "overlap"){
	for ($i = 0; $i <= $alen; $i++){ $DP[$i][0] = 0 }
	for ($j = 0; $j <= $blen; $j++){ $DP[0][$j] = 0 }
    }
    else {
	die;
    }
    return @DP;    
}

sub align_matrix {
    my ($a, $b, $type) = @_;
    my ($i, $j, @DP, @a, $alen, @b, $blen);

    @a = split(//, $a);
    @b = split(//, $b);
    $alen = length($a);
    $blen = length($b);

    @DP = init_matrix($alen, $blen, $type);

    for (my $i = 1; $i <= $alen; $i++){
	for (my $j = 1; $j <= $blen; $j++){
	    my $mat = match($a[$i-1], $b[$j-1]) ? 1 : -1;
	    $DP[$i][$j] = max($DP[$i-1][$j-1] + $mat,
			      $DP[$i-1][$j] + $INDEL, $DP[$i][$j-1] + $INDEL);
	}
    }
    return @DP;
}

sub backtrace_start {
    my ($DP_ref, $type) = @_;
    my ($i, $j, $maxi, $maxj, $max);

    if ($type eq "global"){
	$maxi = scalar(@{$DP_ref})-1;
	$maxj = scalar(@{$DP_ref->[0]}) - 1;
	$max = $DP_ref->[$maxi][$maxj];
    }
    elsif ($type eq "overlap") {
	
	$j = scalar(@{$DP_ref->[0]}) - 1;
	for ($i = 0; $i < scalar(@{$DP_ref}); $i++){
	    if (!defined $max || $max < $DP_ref->[$i][$j]){
		$max = $DP_ref->[$i][$j];
		$maxi = $i;
		$maxj = $j;
	    }
	}
	$i = scalar(@{$DP_ref})-1;
	for ($j = 0; $j < scalar(@{$DP_ref->[0]}); $j++){
	    if (!defined $max || $max < $DP_ref->[$i][$j]){
		$max = $DP_ref->[$i][$j];
		$maxi = $i;
		$maxj = $j;
	    }
	}
    }
    else {
	die;
    }
    return ($maxi, $maxj, $max);
}

sub backtrace_matrix {
    my ($seq0, $seq1, $type, @DP) = @_;
    my ($aln0, $aln1, @seq0, @seq1);

    @seq0 = split(//, $seq0);
    @seq1 = split(//, $seq1);

    $aln0 = $aln1 = "";

    my ($i, $j, $scr) = backtrace_start(\@DP, $type);
    #print $i."\t".$j."\t".$scr."\n";
    if ($i < length($seq0)){
	$aln0 .= reverse(substr($seq0, $i));
	$aln1 .= "-" x (length($seq0) - $i);
    }
    if ($j < length($seq1)){
	$aln0 .= "-" x (length($seq1) - $j);
	$aln1 .= reverse(substr($seq1, $j));
    }
    while($type eq "global" ? 
	  !($i == 0 && $j == 0) : 
	  ($i > 0 && $j > 0)){
	my $mat = match($seq0[$i-1], $seq1[$j-1]) ? 1 : -1 if $i > 0 && $j > 0;
	if ($i > 0 && $j > 0 && $DP[$i][$j] == $DP[$i-1][$j-1] + $mat){
	    $aln0 .= $seq0[--$i];
	    $aln1 .= $seq1[--$j];
	}
	elsif ($i > 0 && $DP[$i][$j] == $DP[$i-1][$j] + $INDEL){
	    $aln0 .= $seq0[--$i];
	    $aln1 .= "-";
	}
	elsif ($j > 0 && $DP[$i][$j] == $DP[$i][$j-1] + $INDEL){
	    $aln0 .= "-";
	    $aln1 .= $seq1[--$j];
	}
	else {
	    die;
	}
	#print $i."\t".$j."\n";
	#print reverse($aln0)."\n";
	#print reverse($aln1)."\n\n";
    }
    if ($i > 0) { 
	$aln0 .= reverse(substr($seq0, 0, $i));
	$aln1 .= "-" x $i;	
    }
    if ($j > 0) {
	$aln0 .= "-" x $j;
	$aln1 .= reverse(substr($seq1, 0, $j));
    }
    $aln0 = reverse($aln0);
    $aln1 = reverse($aln1);
    
    return($scr, $aln0, $aln1);
}

sub match {
    my ($cha, $chb) = @_;
    return 1 if $cha eq $chb;## || $cha eq "N" || $chb eq "N";
    return 0 if $cha ne $chb;
}

sub printAlign {
    my ($aln0, $aln1, $type) = @_;
    my ($inter, @aln0, @aln1);
    die $aln0."\t".$aln1."\n" if length($aln0) != length($aln1);
    @aln0 = split(//, $aln0);
    @aln1 = split(//, $aln1);
    for (my $i = 0; $i < length($aln0); $i++){
	if (match($aln0[$i], $aln1[$i]) &&
	    $aln0[$i] ne "-" && $aln1[$i] ne "-"){
	    $inter .= "|";
	}
	else {
	    $inter .= " ";
	}
    }
    print $aln0."\n".$inter."\n".$aln1."\nedist=".getEdist($aln0, $aln1, $type)."\n";
}

sub getEdist {
    my ($aln0, $aln1, $type) = @_;
    my ($edist, @aln0, @aln1);
    die if length($aln0) != length($aln1);
    if ($type eq "overlap"){
	if ($aln0 =~ m/^([-]+)/ ||
	    $aln1 =~ m/^([-]+)/){
	    $aln0 = substr($aln0, length($1));
	    $aln1 = substr($aln1, length($1));
	}
	if ($aln0 =~ m/([-]+)$/ ||
	    $aln1 =~ m/([-]+)$/){
	    $aln0 = substr($aln0, 0, length($aln0) - length($1));
	    $aln1 = substr($aln1, 0, length($aln1) - length($1));	    
	}
    }

    $edist = 0;
    @aln0 = split(//, $aln0);
    @aln1 = split(//, $aln1);
    for (my $i = 0; $i < length($aln0); $i++){
	if (!match($aln0[$i], $aln1[$i]) ||
	    $aln0[$i] eq "-"){
	    $edist++;
	}
    }
    return $edist;    
}

sub test_align {
    my $a = "ACAGCACTA";
    my $b = "GCA";
    #my $b = "AGCTTCACA";
    my $type = "global";
    my @DP = align_matrix($a, $b, $type);
    #for (my $i = 0; $i < scalar(@DP); $i++){
	#for (my $j = 0; $j < scalar(@{$DP[0]}); $j++){
	#    print $DP[$i][$j]."\t";
	#}
	#print "\n";
    #}
    my ($scr, $aln0, $aln1) = backtrace_matrix($a, $b, $type, @DP);

    print $aln0."\n".$aln1."\nscr=".$scr."\n";
    die;
}

main();
