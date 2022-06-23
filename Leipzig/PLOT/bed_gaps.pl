#!/usr/bin/perl
#perl bed_gaps.pl infile.bed >outfile
# tHis program only cares about the start of the chromosome missing, ends not considered

$file = $ARGV[0];
open(Input,"<","$file");
open (OUT,">","$file.gaps");

@array=<Input>;
$k=0;
print "#query\tstart_gap\tend_gap\tleng_gap\n";

while($array[$k])
{
chomp $array[$k];
chomp $array[$k+1];
@e=split /\t/,$array[$k];
@en=split /\t/,$array[$k+1];
	if($en[0] ne $e[0] && $en[1]>=100)
	{
	$gap=$en[1]-0;
	print OUT "$en[0]\t0\t$en[1]\t$gap\n";
	}

	if($e[0] eq $en[0] && $en[1] > $e[2])
	{
	$gap=$en[1]-$e[2];
	print OUT "$e[0]\t$e[2]\t$en[1]\t$gap\n";
	}

$k++;
}

