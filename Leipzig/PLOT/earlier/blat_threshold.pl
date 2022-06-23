#!/usr/bin/perl

# perl blast_threshold.pl infile > outfile

$filename = $ARGV[0];
$filename_snp = $filename.'_snp';
$filename_per_orth = $filename.'_perfect_orthologs';
$filename_rear = $filename.'_rear';
$filename_unknown = $filename.'_unknown';

open($filehandle, '<', $filename) or die "Could not open $filename\n";
@arr=<$filehandle>;
chomp(@arr);
close(F);
local $" = "\t";			#for printing tab in array


open (fpo, '>>', $filename_per_orth);
open (fsnp, '>>', $filename_snp);
open (frear, '>>', $filename_rear);
open (funkw, '>>', $filename_unknown);


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
print fsnp "@blast\n";
}



#####Definite rearrangements 
elsif ($blast[4] > 0 || $blast[6] > 0)
{
print fsnp "@blast\n";
}



##### Remaining
else
{
	print fremain "@blast\n";
}

}

	close fpo;
	close fsnp;
	close frear;
	close funkw;


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


