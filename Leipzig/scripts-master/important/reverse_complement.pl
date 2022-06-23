#!/usr/bin/env perl

############     reverse_complement.pl     ############
#
#	This Perl script manipulates nucleotide sequences and makes
#	them:
#	- reverse
#	- complementary
#	- reverse + complementary
#	
#	The script is part of the "NGS tools for the novice"
#	authored by David Rosenkranz, Institue of Anthropology
#	Johannes Gutenberg University Mainz, Germany
#
#	Contact: rosenkrd@uni-mainz.de



############     HOW TO USE     ############
#
#	You can pass file names of the files to be processed and
#	the output file name as arguments to the script.
#
#	Input files have to be passed to the scripts via -i:
#	-i file1.fas -i file2.fas
#
#	A bare digit from 1-3 (1 or 2 or 3) will tell the script, how
#	to manipulate your sequence data set:
#	-1 = reverse each sequence.
#	-2 = make each sequence complementary.
#	-3 = make each each sequence reverse complementary.
#	By default, the value is set to 3.
#
#	For example you can type the following command:
#	perl reverse_complement.pl -i file1.fas -3
#
#	Output files will have the name of the input file extended
#	by _r, _c or _rc.
#	If your input file is input.fas, the output files will be
#	input_r.fas
#	input_c.fas
#	input_rc.fas
#
#	If you do not want to enter each file name seperately, you
#	can provide a file that contains a list all file names (one
#	file name per line). Pass the file name to the script via -I:
#	-I list_of_files.txt
#
#	Multiple files and combinations of all kinds of arguments are
#	allowed:
#	perl reverse_complement.pl -i file1.fas -I list_of_files.txt -3




@input_files=();
$manipulate=3;
$|=1;

###   CHECK COMMAND LINE ARGUMENTS   ###
if(@ARGV==0)
	{
	print"No arguments passed to the script!\nIf you entered arguments try the following command:\nperl reverse_complement.pl -argument1 -argument2 ...\n\n";
	exit;
	}

$argv="";
foreach(@ARGV)
	{
	$argv.=$_;
	}
@arguments=split('-',$argv);

foreach(@arguments)
	{
	if($_=~/^ *i/)
		{
		$_=~s/^ *i//;
		$_=~s/ //g;
		push(@input_files,$_);
		}
	elsif($_=~/^ *I/)
		{
		$_=~s/^ *I//;
		$_=~s/ //g;
		open(FILES_IN,"$_");
		while(<FILES_IN>)
			{
			unless($_=~/^\s*$/)
				{
				$_=~s/\s//sg;
				push(@input_files,$_);
				}
			}
		}
	elsif($_=~/^ *[123] *$/)
		{
		$_=~s/ *//g;
		$_=~s/ //g;
		$manipulate=$_;
		}
	elsif($_!~/^\s*$/)
		{
		print"Don't know how to treat argument $_!\nIt will be ignored.\n\n";
		}
	}
if(@input_files==0)
	{
	print"No input file specified!\n";
	exit;
	}

###   PRINT ARGUMENTS   ###
print"The following files will be processed:\n";
foreach(@input_files)
	{
	if(-e $_)
		{
		print"$_\n";
		push(@input_files_ok,$_);
		}
	else
		{
		print"could not find file: $_. It will be ignored.\n";
		}
	}
print"\nSequences will be made: ";
if($manipulate==1)
	{
	print"REVERSE.\n\n";
	}
elsif($manipulate==2)
	{
	print"COMPLEMENTARY.\n\n";
	}
elsif($manipulate==3)
	{
	print"REVERSE COMPLEMENTARY.\n\n";
	}

###   START   ###
$number_of_outfiles=0;
foreach$file(@input_files_ok)
	{
	$number_of_outfiles++;
	$out_file_name=$file;
	$out_file_name=~s/\..+$//;
	$out_file_name.="($number_of_outfiles).fas";
	open(OUT,">$out_file_name");
	print"processing $file";
	open(IN,$file);
	while(<IN>)
		{
		if($_=~/^>/)
			{
			$title=$_;
			}
		elsif($_!~/^\s*$/)
			{
			$seq=$_;
			chomp$seq;
			if($manipulate==1)
				{
				$seq=reverse$seq;
				}
			if($manipulate==2)
				{
				$complementary=$seq;
				$complementary=~tr/ATGCatgc/TACGtacg/;
				$seq=$complementary;
				}
			if($manipulate==3)
				{
				$reverse_complementary=reverse$seq;
				$reverse_complementary=~tr/ATGCatgc/TACGtacg/;
				$seq=$reverse_complementary;
				}
			print OUT "$title$seq\n";
			}
		}
	close IN;
	close OUT;
	print" done.\n" 
	}
exit;