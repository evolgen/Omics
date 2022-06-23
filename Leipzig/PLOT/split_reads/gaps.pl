#!/usr/bin/perl
$file = $ARGV[0];
open(Input,"<","$file");
open (OUT,">","$file.gaps");
@array=<Input>;
$k=0;
print OUT "#query\tstart_gap\tend_gap\tleng_gap\tsign\n";
while($array[$k])
{
chomp $array[$k];
chomp $array[$k+1];
@e=split /\t/,$array[$k];
@en=split /\t/,$array[$k+1];
	if($e[0] eq $en[0] && $e[3] eq $en[3])
	{
	$gap=$en[1]-$e[2];
	print OUT "$e[0]\t$e[2]\t$en[1]\t$gap\t$e[3]\n";
	}
$k++;
}

