#!/usr/bin/env perl
$file = $ARGV[0];
$fileout = $ARGV[1];
open(Input,"<","$file");
open (OUT,">","$fileout");
@array=<Input>;
$k=0;$r=1;

while($array[$k] && $array[$r])
{
	chomp $array[$k];
	chomp $array[$r];
	@e=split /\t/,$array[$k];
	@en=split /\t/,$array[$r];

if($e[0] eq $en[0] && $e[2] ge $en[1] && $e[6] eq $en[6]) {
	print OUT "$array[$k]\n$array[$r]\n";
	$r++;$k++;
}

elsif($en[1]>=$e[1] &&  $en[2]<=$e[2] && $en[6]<=$e[6]) {
	$r++;
}

else {
$k++;$r++;
}

}

