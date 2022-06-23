#!/usr/bin/perl

# perl blast_threshold.pl infile > outfile

$filename = $ARGV[0];

$filename_snp = $filename.'_snp';
$filename_per_orth = $filename.'_perfect_orthologs';
$filename_single_idt_5bp = $filename.'_single_cr_5bp';
$filename_single_idt_16bp = $filename.'_single_cr_16bp';
$filename_single_idt_30bp = $filename.'_single_cr_30bp';
$filename_single_idt_large = $filename.'_single_cr_large';
$filename_single_idt_uneven = $filename.'_single_cr_uneven';

$filename_multiple_idt_inone = $filename.'_multiple_cr_inone';
$filename_single_cr_both = $filename.'_single_cr_both';
$filename_multi_exact_cr = $filename.'_multiple_exact_cr';
$filename_multi_uneven_cr = $filename.'_multiple_uneven_cr';

#$filename_single_inv = $filename.'_single_inv';
$filename_remaining = $filename.'_remaining';

open($filehandle, '<', $filename) or die "Could not open $filename\n";
@arr=<$filehandle>;
chomp(@arr);
close(F);
local $" = "\t";			#for printing tab in array

open (fpo, '>>', $filename_per_orth);
open (fsidt5, '>>', $filename_single_idt_5bp);
open (fsidt16, '>>', $filename_single_idt_16bp);
open (fsidt30, '>>', $filename_single_idt_30bp);
open (fsidtL, '>>', $filename_single_idt_large);
open (fsidtue, '>>', $filename_single_idt_uneven);
open (fsnp, '>>', $filename_snpdata);
open (fmidtio, '>>', $filename_multiple_idt_inone);
open (fscrb, '>>', $filename_single_cr_both);
open (fmecr, '>>', $filename_multi_exact_cr);
open (fmuecr, '>>', $filename_multi_uneven_cr);
open (fremain, '>>', $filename_remaining);

for ($j=0;$j<=$#arr;$j++)
{
@blast=split(/\t/,$arr[$j]);

#####Perfect orthologs
if ($blast[4] == 0 && $blast[6] == 0  &&  $blast[1] == 0)
{
print fpo "@blast\n";
}


#####Definite Polymorphism and SNP
elsif ($blast[4] == 0 && $blast[6] == 0 && $blast[1] > 0)
{
#open (fsnp, '>>', $filename_snpdata);
print fsnp "@blast\n";
}


######Single CR in both
elsif ($blast[4] == 1 && $blast[6] == 1)
{
#open (fscrb, '>>', $filename_single_cr_both);
print fscrb "@blast\n";
}


#####Definite single CR in one, based on sizes
elsif (($blast[4] == 1 && $blast[6] == 0) || ($blast[4] == 0 && $blast[6] == 1))				#Only one gap i.e. either in query or target
{
	if ($blast[5] <= 5 && $blast[7] <= 5)
	{
#open (fsidt5, '>>', $filename_single_idt_5bp);
	print fsidt5 "@blast\n";
	}
#	close fsidt5;
	elsif (($blast[5] > 5 && $blast[5] <= 16) && ($blast[7] > 5 && $blast[7] <= 16))
	{  
#open (fsidt16, '>>', $filename_single_idt_16bp);
	print fsidt16 "@blast\n";
	}
#	close fsidt16;
	elsif (($blast[5] > 16 && $blast[5] <= 30) && ($blast[7] > 16 && $blast[7] <= 30))
	{  
#open (fsidt30, '>>', $filename_single_idt_30bp);
	print fsidt30 "@blast\n";
	}
#	close fsidt30;

	elsif (($blast[5] > 30) && ($blast[7] > 30))
	{  
#open (fsidtL, '>>', $filename_single_idt_large);
	print fsidtL "@blast\n";
	}
#	close fsidtL;
	else
	{  
#open (fsidtL, '>>', $filename_single_idt_large);
	print fsidtue "@blast\n";
	}
#	close fsidtL;
}


#####Multiple indels-transpositions
elsif (($blast[4] > 1 && $blast[6] == 0) || ($blast[4] == 0 && $blast[6] > 1))
{
#open (fmidtio, '>>', $filename_multiple_idt_inone);
print fmidtio "@blast\n";
}


#####Multiple CRs in both
elsif (($blast[4] > 1 && $blast[6] > 1))
{
	if ($blast[4] = $blast[6])
	{
#open (fmecr, '>>', $filename_multi_exact_cr);
	print fmecr "@blast\n";
	}
#	close fmecr;
	if ($blast[4] != $blast[6])
	{
#open (fmuecr, '>>', $filename_multi_uneven_cr);
	print fmuecr "@blast\n";
	}
#	close fmuecr;
}


##### Remaining results
else
{
#open (fremain, '>>', $filename_remaining);
print fremain "@blast\n";
}

}
	close fsidt5;
	close fsidt16;
	close fsidt30;
	close fsidtL;
	close fsidtue;
	close fsnp;
	close fmidtio;
	close fscrb;
	close fmecr;
	close fmuecr;
	close fremain;


=begin comment
1. match           number of matching nucleotides in alignment that aren't part of repeats
2. mismatch        number of nucleotides in alignment that don't match
3. rep. match      number of nucleotides in alignment that are part of repeats
4. N's             number of N's in alignment
5. Q gap count     number of inserts in query sequence
6. Q gap bases     number of nucleotides inserted in query sequence
7. T gap count     number of inserts in target sequence (chromosome)
8. T gap bases     number of nucleotides inserted in target sequence (chromosome)
9. strand          chromsome strand (+ or -)
10. Q name          name of query sequence
11. Q size          length of query sequence
12. Q start         start of query sequence in alignment
13. Q end           end of query sequence in alignment
14. T name          matching chromsome ("target", ex: chr13)
15. T size          length of target sequence (chromosome)
16. T start         start of target sequence (chromosome) comprising alignment
17. T end           end of target sequence (chromosome) comprising alignment
18. block count     number of blocks (may be exons for cDNA) of matching regions
19. blockSizes      sizes of blocks (may be exons for cDNA) of matching regions (delimited by commas)
20. qStarts         list of query nts at starts of blocks (delimited by commas)
21. tStarts         list of chromsome nts at starts of blocks (delimited by commas)
=end comment


