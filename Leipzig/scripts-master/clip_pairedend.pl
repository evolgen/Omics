#!/usr/bin/perl -w
use strict;
use File::Basename;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Getopt::Std;
use List::Util;
use Cwd;
use IO::File;
use POSIX qw(tmpnam);
use FindBin qw($Bin);
use threads;


use vars qw ($help $mate1 $mate2 $out1 $out2 $adapter1 $adapter2 $Ns $stat1 $stat2);
#my $cutadapt = "/scr/k41san/tools/cutadapt/cutadapt-1.3/bin/cutadapt";
my $cutadapt = "cutadapt";
$Ns = 1;

# -----------------------------------------------------------------------------
# OPTIONS

GetOptions (
"m1=s"       => \$mate1,
"m2=s"       => \$mate2,
"o1=s"       => \$out1,
"o2=s"       => \$out2,
"a1=s"       => \$adapter1,
"a2=s"       => \$adapter2,
"s1=s"       => \$stat1,
"s2=s"       => \$stat2,
"n=s"        => \$Ns,
"help"       => \$help,
"h"          => \$help);
usage() if ($help || !$mate1 || !$mate2 || !$out1 || !$out2 || !$stat1 || !$stat2);

# MAIN

my ($name1, $type1, $gz1) = getName($mate1);
my ($name2, $type2, $gz2) = getName($mate2);

my $thread1 = threads->create('doCall', "$cutadapt -e 0.15 -m 0 -a $adapter1 -o $name1.tmp$type1$gz1 $mate1 > $stat1");
my $thread2 = threads->create('doCall', "$cutadapt -e 0.15 -m 0 -a $adapter2 -o $name2.tmp$type2$gz2 $mate2 > $stat2");

$thread1->join();
$thread2->join();

if($gz1 eq ".gz"){
    open(MATE1, "zcat $mate1 | paste - - - - | ") || die "cannot open $mate1\n";
    open(CLIPPED1, "zcat $name1.tmp$type1$gz1 | paste - - - - | ") || die "cannot open $name1.tmp$type1$gz1\n";
}else{
    open(MATE1, "cat $mate1 | paste - - - - | ") || die "cannot open $mate1\n";
    open(CLIPPED1, "cat $name1.tmp$type1$gz1 | paste - - - - | ") || die "cannot open $name1.tmp$type1$gz1\n";
}
if($gz2 eq ".gz"){
    open(MATE2, "zcat $mate2 | paste - - - - | ") || die "cannot open $mate2\n";
    open(CLIPPED2, "zcat $name2.tmp$type2$gz2 | paste - - - - | ") || die "cannot open $name2.tmp$type2$gz2\n";
}else{
    open(MATE2, "cat $mate2 | paste - - - - | ") || die "cannot open $mate2\n";
    open(CLIPPED2, "cat $name2.tmp$type2$gz2 | paste - - - - | ") || die "cannot open $name2.tmp$type2$gz2\n";
}

if($out1 !~ /\.gz/){$out1 .= ".gz";}
if($out2 !~ /\.gz/){$out2 .= ".gz";}
open(OUT1, " | gzip -c > $out1") || die "cannot open $out1\n";
open(OUT2, " | gzip -c > $out2") || die "cannot open $out2\n";
while (<MATE1>){
    my $m1 = $_;
    chomp($m1);
    my @m1e = split(/\t/,$m1);
    my $m2 = <MATE2>;
    chomp($m2);
    my @m2e = split(/\t/,$m2);
    my $c1 = <CLIPPED1>;
    chomp($c1);
    my @c1e = split(/\t/,$c1);
    my $c2 = <CLIPPED2>;
    chomp($c2);
    my @c2e = split(/\t/,$c2);
    
    # both mates are too short
    next if(length($c1e[1]) < 15 && length($c2e[1]) < 15);
    
    # mate1 is too short
    if(length($c1e[1]) < 15 && $Ns == 0){
        print OUT1 join("\n", @m1e)."\n";
        print OUT2 join("\n", @c2e)."\n";
        next;
    }
    if(length($c1e[1]) < 15 && $Ns == 1){
        $c1e[1] = "N"; $c1e[3] = "!";
        print OUT1 join("\n", @c1e)."\n";
        print OUT2 join("\n", @c2e)."\n";
        next;
    }
    
    # mate2 is too short
    if(length($c2e[1]) < 15 && $Ns == 0){
        print OUT1 join("\n", @c1e)."\n";
        print OUT2 join("\n", @m2e)."\n";
        next;
    }
    if(length($c2e[1]) < 15 && $Ns == 1){
        $c2e[1] = "N"; $c2e[3] = "!";
        print OUT1 join("\n", @c1e)."\n";
        print OUT2 join("\n", @c2e)."\n";
        next;
    }
    
    #both mates are correct
    print OUT1 join("\n", @c1e)."\n";
    print OUT2 join("\n", @c2e)."\n";
}
close(MATE1);
close(MATE2);
close(CLIPPED1);
close(CLIPPED2);
close(OUT1);
close(OUT2);

system("rm $name1.tmp$type1$gz1");
system("rm $name2.tmp$type2$gz2");

# -----------------------------------------------------------------------------
sub usage {
    print STDERR "\nusage: clipPairedEnd.pl -m1 <file> -m2 <file> -o1 <file> -o2 <file> -s1 <file> -s2 <file>\n";
    print STDERR "overlap reads with annotations\n";
    print STDERR "\n";
    print STDERR "[INPUT]\n";
    print STDERR " -m1 <file>    raw mates 1\n";
    print STDERR " -m2 <file>    raw mates 2\n";
    print STDERR " -a1 <string>  adapter for mates 1\n";
    print STDERR " -a2 <string>  adapter for mates 2\n";
    print STDERR " -o1 <file>    clipped mates 1\n";
    print STDERR " -o2 <file>    clipped mates 2\n";
    print STDERR " -s1 <file>    clippStat mates 1\n";
    print STDERR " -s2 <file>    clippStat mates 2\n";
    print STDERR " -n  <int>     1: fill mates <15nt with Ns (default)\n";
    print STDERR "               0: reset mates <15nt with original mate\n";
    print STDERR " -h <file>     this (usefull) help message\n";
    print STDERR "[VERSION]\n";
    print STDERR " 12-21-2012\n";
    print STDERR "[BUGS]\n";
    print STDERR " Please report bugs to david\@bioinf.uni-leipzig.de\n";
    print STDERR "\n";
    exit(-1);
}

sub getName{
    my ($file) = @_;
    my $gz = ""; my $type = ".fq";
    my @folderStruct = split(/\//,$file);
    my $name = $folderStruct[$#folderStruct];
    if($name =~ /\.gz/){
        $name =~ s/\.gz//;
        $gz = ".gz";
    }
    if($name =~ /\.fa$/){
        $name =~ s/\.fa//;
        $type = ".fa";
    }
    $name =~ s/\.fq//;
    $name =~ s/\.fastq//;
    return ($name, $type, $gz);
}

sub doCall{
    # make the calling
    my ($string) = @_;
    system("$string");
}
