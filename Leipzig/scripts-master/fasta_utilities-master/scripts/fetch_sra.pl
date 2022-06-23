#!/usr/bin/perl
use warnings;
use strict;
###############################################################################
# By Jim Hester
# Date:05/15/2012
# Last Modified: 2013 Jul 11 02:46:22 PM
# Title:fetch_sra.pl
# Purpose:this script downloads the sra sequences from NCBI using aspera and outputs a fastq file
###############################################################################

###############################################################################
# Code to handle help menu and man page
###############################################################################
use Getopt::Long;
use Pod::Usage;
my $man = 0;
my $help = 0;
my $asperaDir = "~/.aspera";
my $fastq_dump_options = "-F --gzip --split-3";
my $fastq_dump_dir = "~/packages/sratoolkit.2.3.2-4-centos_linux64/bin";
GetOptions('aspera=s' => \$asperaDir, 'fastq-dump=s' => \$fastq_dump_dir,
'help|?' => \$help, man => \$man) or pod2usage(2); pod2usage(2) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage("$0: No files given.")  if ((@ARGV == 0) && (-t STDIN));

######################c########################################################
# Automatically extract compressed files
###############################################################################
@ARGV = map { s/(.*\.gz)\s*$/pigz -dc < $1|/; s/(.*\.bz2)\s*$/pbzip2 -dc < $1|/;$_ } @ARGV;
###############################################################################
# fetch_sra.pl
###############################################################################

my $link = shift;
my $description = shift;

if($link =~ m{ftp://([^/]+)(/.+/)([^/]+)}){
  my($server,$path,$dataset)=($1,$2,$3);
  my $command = "$asperaDir/connect/bin/ascp -o Overwrite=never -k2 -i $asperaDir/connect/etc/asperaweb_id_dsa.putty -QT -L . -l 500m anonftp\@$server:$path$dataset .";
  print STDERR $command, "\n";
  if(not -e $dataset){
    system($command);
    system("ln -s $dataset $description") if $description;
  }
  for my $file (glob "$dataset/*sra*"){
    my $dump_command = "$fastq_dump_dir/fastq-dump $fastq_dump_options $file";
    print STDERR $dump_command, "\n";
    system($dump_command) unless glob "$dataset*fastq*";
  }
}
###############################################################################
# Help Documentation
###############################################################################

=head1 NAME

fetch_sra.pl - this script downloads the sra sequences from NCBI using aspera and outputs a fastq file

=head1 SYNOPSIS

fetch_sra.pl [options] [file ...]

Options:
      -help
      -man               for more info

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<fetch_sra.pl> this script downloads the sra sequences from NCBI using aspera and outputs a fastq file

=cut

