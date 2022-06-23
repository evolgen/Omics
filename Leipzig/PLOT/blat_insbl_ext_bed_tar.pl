#!/usr/bin/perl 
### Script to extract the insertion gaps only

use strict;
my $inputfile=$ARGV[0];
if (! -e $inputfile) {die "Please  provide an input blat file\n";}

my $outfile = $inputfile;
$outfile =~ s/\.psl*//;
my $outfile1 = $outfile."_query_gap.bed";
my $outfile2 = $outfile."_target_gap.bed";

#system ("rm $outfile1_falsegaps");


if ( -e $outfile2)      {
system("rm $outfile2");
}

#if ( -e $outfile1 )      {
#system("rm $outfile1");
#}



open (my $info, "<", $inputfile) or die "Could not open $inputfile\n";
##open (my $outQuery, ">", $outfile1) or die "Could not create $outfile1\n";  
open (my $outTarget, ">", $outfile2) or die "Could not create $outfile2\n";


my $linenumber = 0;
my $blockis = 1;
while( my $theline = <$info>)  {   
	chomp $theline;
	$linenumber = $linenumber+1;
	my @thisline = split("\t", $theline);
	if (scalar(@thisline) != 21) {die "Please provide a valid blat output file in .psl format\n"};
	my $refname = $thisline[9];	
	my $blockcount = $thisline[17];
	my @blocksizes = split(",", $thisline[18]);
	my @qStarts = split(",", $thisline[19]);
	my @tStarts = split(",", $thisline[20]);
	my $Line_count = 1;

	for (my $count=0; $count < $blockcount-1; $count++) {
		my $blocknumber = $count+1;
		if ($thisline[8] eq '+') {
#		if ($thisline[8] eq '+')	{
#10,11,12			14,15,16
			my $count2=$count+1;
                        my $QueryEnd = $qStarts[$count2];
			my $QueryStart = $qStarts[$count] + $blocksizes[$count];
			my $InsertsizeQ = $QueryEnd - $QueryStart;
                        my $TargetEnd = $tStarts[$count2];
                        my $TargetStart = $tStarts[$count] + $blocksizes[$count];
			my $InsertsizeT = $TargetEnd - $TargetStart;
			##print $outQuery $thisline[9],"\t",$QueryStart,"\t",$QueryEnd,"\t",$thisline[10],"\t",$thisline[13],"\t+\t",$thisline[10],"+",$blocknumber,":",$InsertsizeQ,"\t",$thisline[13],":",$InsertsizeT,"\t",$InsertsizeQ,"\t",$InsertsizeT,"\n";
#			print $outTarget $thisline[13],"\t",$TargetStart,"\t",$TargetEnd,"\t",$thisline[14],"\t",$thisline[14],"+",$blocknumber,":",$InsertsizeT,"\t+\t",$thisline[9],"\t",$thisline[10],"\t",$thisline[9],"+",$blocknumber,":",$InsertsizeQ,"\t",$TargetStart,"\t",$TargetEnd,"\t",$InsertsizeT,"\t",$QueryStart,"\t",$QueryEnd,"\t",$InsertsizeQ,"\n";
			print $outTarget $thisline[13],"\t",$TargetStart,"\t",$TargetEnd,"\t",$thisline[14],"\t",$thisline[13],"+",$blocknumber,":",$InsertsizeT,"\t+\t",$thisline[9],"\t",$QueryStart,"\t",$QueryEnd,"\t",$thisline[10],"\t",$linenumber,"_#_",$Line_count,"\t",$InsertsizeT,"\t",$InsertsizeQ,"\n";
			}
# O_15k_250       13306   13381   16556   16556+1:75      +       O_15k_250_@del  16481   O_15k_250_@del+1:0      13306   13381   75      13306   13306   0


		elsif ($thisline[8] eq '-') {
                        my $count2=$count+1;
                        my $QueryEnd = $qStarts[$count2] - 1;
			my $QueryStart = $qStarts[$count] + $blocksizes[$count];
			my $InsertsizeQ = $QueryEnd - $QueryStart + 1;
                        my $TargetEnd = $tStarts[$count2];
                        my $TargetStart = $tStarts[$count] + $blocksizes[$count];
			my $InsertsizeT = $TargetEnd - $TargetStart;
			##print $outQuery $thisline[9],"\t",$QueryStart,"\t",$QueryEnd,"\t",$thisline[10],"\t",$thisline[13],"\t-\t",$thisline[10],"-",$blocknumber,":",$InsertsizeQ,"\t",$thisline[13],":",$InsertsizeT,"\t",$InsertsizeQ,"\t",$InsertsizeT,"\n";
#			print $outTarget $thisline[13],"\t",$TargetStart,"\t",$TargetEnd,"\t",$thisline[14],"\t",$thisline[14],"+",$blocknumber,":",$InsertsizeT,"\t-\t",$thisline[9],"\t",$thisline[10],"\t",$thisline[9],"-",$blocknumber,":",$InsertsizeQ,"\t",$TargetStart,"\t",$TargetEnd,"\t",$InsertsizeT,"\t",$QueryStart,"\t",$QueryEnd,"\t",$InsertsizeQ,"\n";
			print $outTarget $thisline[13],"\t",$TargetStart,"\t",$TargetEnd,"\t",$thisline[14],"\t",$thisline[13],"+",$blocknumber,":",$InsertsizeT,"\t-\t",$thisline[9],"\t",$QueryStart,"\t",$QueryEnd,"\t",$thisline[10],"\t",$linenumber,"_#_",$Line_count,"\t",$InsertsizeT,"\t",$InsertsizeQ,"\n";
				}

		else	{
			print "Error in Strand in line $linenumber as Strand is $thisline[8]\n";
			exit;
			}

		$Line_count++;
		}
	}


close $info;
#close $outQuery;
close $outTarget;


=head
# bedtools getfasta -s -fi querysequence.fa -name -bed query_bedfile.bed -s -fo query_insert_extract.fa
# bedtools getfasta -s -fi targetsequence.fa -name -bed target_bedfile.bed -s -fo target_insert_extract.fa

155	9	0	0	3	18	2	21	+	clint_64075375	200	18	200	vailant_86051643	258	13	198	4	27,7,75,55,	18,53,61,145,	13,55,62,143,

0	matches - Number of bases that match that aren't repeats
1	misMatches - Number of bases that don't match
2	repMatches - Number of bases that match but are part of repeats
3	nCount - Number of 'N' bases
4	qNumInsert - Number of inserts in query
5	qBaseInsert - Number of bases inserted in query
6	tNumInsert - Number of inserts in target
7	BaseInsert - Number of bases inserted in target
8	strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
9	qName - Query sequence name
10	qSize - Query sequence size
11	qStart - Alignment start position in query
12	qEnd - Alignment end position in query
13	tName - Target sequence name
14	tSize - Target sequence size
15	tStart - Alignment start position in target
16	tEnd - Alignment end position in target
17	blockCount - Number of blocks in the alignment (a block contains no gaps)
18	blockSizes - Comma-separated list of sizes of each block
19	qStarts - Comma-separated list of starting positions of each block in query
20	tStarts - Comma-separated list of starting positions of each block in target

