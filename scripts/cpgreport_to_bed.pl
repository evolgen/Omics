#!/usr/bin/env perl

#use warnings;
#use strict;
use Data::Dumper qw(Dumper);

#my @sequence_name;
#my $num1;
#my %hash1;

#my $num_args = $#ARGV + 1;
#if ($num_args != 1) {
if ($#ARGV != 0 ) {
    print "\nUsage: perl ~/RGP/scripts/cpgreport_to_bed.pl $file1\n";
    exit;
}

my $file1 = $ARGV[0];
chomp($file1); 

open(my $infohd1, '<:encoding(UTF-8)', $file1)
    or die "Could not open $file1: $!";

my $current_ID ; 
while( my $line1 = <$infohd1>)  {   
    chomp($line1);
    if($line1 =~ /^ID::/) {
        $current_ID = $line1;
        $current_ID =~ s/^ID:://;
    }
    elsif($line1 =~ /^pos::/) {
        $line1 =~ s/^pos:://;
        print "$current_ID\t$line1\t+\n";
    }
}


close $info;



