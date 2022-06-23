#!/usr/bin/env perl
#perl raplace_with_hash.pl patternfile infile > outfile
use warnings;
use strict;

my $hashfile = $ARGV[0];
my $matchfile = $ARGV[1];

open (my $EQ, '<', $hashfile) or die "cannot find $hashfile";
my %subst;
while (<$EQ>) {
    chomp;                                   
    my ($search, $replace) = split /=/;
    $subst{$search} = $replace;
}

my $regex = join '|', map quotemeta, keys %subst;

open (my $LST, '<', $matchfile) or die "cannot find $matchfile";
while (<$LST>) {
    s/($regex)/$subst{$1}/g;
    print;
}
