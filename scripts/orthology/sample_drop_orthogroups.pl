#!/usr/bin/env perl
use strict;
use warnings;

#Usage:  ls | grep FNA  | perl /global/scratch/gregoryowens/sebastes/bin/find_sample_to_drop.pl > dropping_list.txt
my @keep_species = qw(Sebastes_aleutianus Sebastes_entomelas Sebastes_miniatus Sebastes_pinniger Sebastes_rosaceus Sebastes_umbrosus);
my %keep_species = map { $_ => 1 } @keep_species;
my %hash;
my %reverse_hash;

#pipe in a list of files;
while(<STDIN>){
  chomp;
  my $file = $_;
  open FILE, $file;
  while(<FILE>){
    chomp;
    if ($_ =~ m/^>/){
      my $s = $_;
      $s =~ s/>//g;
      my $sample = $s;
      $hash{$sample}{$file}++;
      $reverse_hash{$file}{$sample}++;
    }
  }
  close FILE;
  print STDERR "Loaded $file\n";
}

my $max_genes = scalar keys %reverse_hash;
print STDERR "Max number of genes is $max_genes\n";
my $max_samples = scalar keys %hash;
my %removed_samples;
my @removed_samples;

#Count the number of genes with full data in all $samples
my $original_kept_genes;
foreach my $gene (sort keys %reverse_hash){
  my $species_sequenced = 0;
  foreach my $iterating_sample (sort keys %{$reverse_hash{$gene}}){
    $species_sequenced++;
  }
  if ($species_sequenced == $max_samples){
    $original_kept_genes++;
  }
}
print STDERR "Iterating over posible samples to remove\n";
foreach my $i (1..($max_samples -6)){
#foreach my $i (1..2){
  my $n_retained = $max_samples - $i;
#  print STDERR "Max number of samples is $n_retained\n";
  #Test to see how many genes would retain if this sample were dropped.
  my %kept_genes;
  foreach my $sample (sort keys %hash){
    if ($keep_species{$sample}){next;}
    if ($removed_samples{$sample}){next;}
    foreach my $gene (sort keys %reverse_hash){
      my $species_sequenced = 0;
      foreach my $iterating_sample (sort keys %{$reverse_hash{$gene}}){
        if ($iterating_sample eq $sample){next;}
        if ($removed_samples{$iterating_sample}){next;}
        $species_sequenced++;
      }
      if ($species_sequenced == $n_retained){
        $kept_genes{$sample}++;
      }
    }
#    print STDERR "If you remove $sample you get $kept_genes{$sample} genes\n";
  }
  #Find sample to drop
  my @ordered_samples = (sort { $kept_genes{$b} <=> $kept_genes{$a} } keys %kept_genes);
  my $dropped_sample = $ordered_samples[0];
  print STDERR "Drop sample $dropped_sample to retain $kept_genes{$dropped_sample} genes\n";
  push (@removed_samples, $dropped_sample);
  $removed_samples{$dropped_sample} = $kept_genes{$dropped_sample};
}
print "round\tsample\tretained_genes";
print "\n0\tNA\t$original_kept_genes";
foreach my $i (0..$#removed_samples){
  my $round = $i+1;
  print "\n$round\t$removed_samples[$i]\t$removed_samples{$removed_samples[$i]}";
}
