#!/usr/bin/perl 

use strict;
my $inputfile=$ARGV[0];
if (! -e $inputfile) {
	die "Please  provide an input blat file\n";	}

#my $outfile = $inputfile;
#$outfile =~ s/\.psl*//;
#my $outfile1 = $outfile."_filt_overlap.bed";

#system ("rm $outfile1_falsegaps");
local $" = "\t";

#my @data = <$inputfile>;
my %set_hash;

#while (my $line = <@data>) {
while (my $line = <$inputfile>) {
    chomp $line;
    my @cols = split(/\t/, $line);

    # if there is already a value with key $cols[0] in the hash
    if (exists $set_hash{$cols[0]}) {
	print "$cols[0] already exists in the hash\n";
    }
    
# otherwise just add the data to the hash
    else {
	$set_hash{$cols[0] = $cols[1];
    }
}
}





#==head
# bedtools getfasta -s -fi querysequence.fa -name -bed query_bedfile.bed -s -fo query_insert_extract.fa
# bedtools getfasta -s -fi targetsequence.fa -name -bed target_bedfile.bed -s -fo target_insert_extract.fa

#155	9	0	0	3	18	2	21	+	clint_64075375	200	18	200	vailant_86051643	258	13	198	4	27,7,75,55,	18,53,61,145,	13,55,62,143,

#0	matches - Number of bases that match that aren't repeats
#1	misMatches - Number of bases that don't match
#2	repMatches - Number of bases that match but are part of repeats
#3	nCount - Number of 'N' bases
#4	qNumInsert - Number of inserts in query
#5	qBaseInsert - Number of bases inserted in query
#6	tNumInsert - Number of inserts in target
#7	BaseInsert - Number of bases inserted in target
#8	strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
#9	qName - Query sequence name
#10	qSize - Query sequence size
#11	qStart - Alignment start position in query
#12	qEnd - Alignment end position in query
#13	tName - Target sequence name
#14	tSize - Target sequence size
#15	tStart - Alignment start position in target
#16	tEnd - Alignment end position in target
#17	blockCount - Number of blocks in the alignment (a block contains no gaps)
#18	blockSizes - Comma-separated list of sizes of each block
#19	qStarts - Comma-separated list of starting positions of each block in query
#20	tStarts - Comma-separated list of starting positions of each block in target
