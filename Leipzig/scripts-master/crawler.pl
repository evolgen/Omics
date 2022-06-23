#!/usr/bin/perl

=head1 DESCRIPTION

	crawler.pl

=head1 USAGE

	perl crawler.pl

=head2 Options


	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message. press q to exit this screen.	

=head1 Author

	Sunit Jain, (Fri Sep  6 13:48:16 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;
use File::Spec;
use File::Path;
use File::Copy;

my $version="crawler.pl\tv0.0.1b";
my $list="contents.list";
my $path= "/geomicro/data1/COMMON/scripts/";
my $sandbox="/geomicro/data1/COMMON/scripts/sandbox";

GetOptions(
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);
print "\# $version\n";

open(LIST, "<".$list) || die $!;
my $folder="./ERROR/";
mkpath($folder);
while(my $line=<LIST>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;
	
	if ($line=~ /\:$/){
	# Folders
		$line=~ s/\:$//;
		$folder=File::Spec->catfile($path, $line);
		if (! -d $folder){
			my $dir = eval{ mkpath($folder) };
			die "Failed to create $line: $@\n" unless $dir;
		}
	}
	else{
	# Files
		my($fileName, @comment)=split(/\t/, $line);
		my $file=File::Spec->catfile($sandbox, $fileName);
		unless($file=~ /\./g){
			warn "'$file' must be a folder. Skipping...\n";
			next;
		}
		
		if(! -e $file){
			warn "Not Found: $file\n";
			next;
		}
		copy($file, $folder) || die "Failed to copy $file: $!\n";
	}
}
close LIST;
