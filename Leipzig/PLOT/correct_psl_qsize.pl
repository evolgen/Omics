#!/usr/bin/env perl
# perl split_hit.pl infile >outfile
# infile pattern
# Query	Start	Stop	Length	Target Strand	Target_len	Start	Stop
# If psl, then <(awk -F'\t' '{print $10"\t"$12"\t"$13"\t"$11"\t"$14"\t"$9"\t"$15"\t"$16"\t"$17}' infile | sort -k1,1 -k6,6 -k2,2n)
#
#
#
$file = $ARGV[0];
$fileout = $ARGV[1];
open(Input,"<","$file");

#open (OUT,">","$fileout");

@array=<Input>;
$k=0;$r=1;

while($array[$k] && $array[$r])
{
	chomp $array[$k];
	chomp $array[$r];
	@e=split /\t/,$array[$k];
	@en=split /\t/,$array[$r];

# First line
if($k==0) {
        print "$array[$k]\n";
	$val=$e[10];
	if ($e[9] eq $en[9]) {
        print "$en[0]\t$en[1]\t$en[2]\t$en[3]\t$en[4]\t$en[5]\t$en[6]\t$en[7]\t$en[8]\t$en[9]\t$val\t$en[11]\t$en[12]\t$en[13]\t$en[14]\t$en[15]\t$en[16]\t$en[17]\t$en[18]\t$en[19]\t$en[20]\t$en[21]\n";
#	print "$e[10]\n";
	}
	else {
		$val=$en[10];
		print "$array[$r]\n";
	}
	$r++;$k++;  
}

# Following lines
# If equal
else {
	if ($e[9] eq $en[9]) {
#	if ($e[10] ne $en[10]) {
        print "$en[0]\t$en[1]\t$en[2]\t$en[3]\t$en[4]\t$en[5]\t$en[6]\t$en[7]\t$en[8]\t$en[9]\t$val\t$en[11]\t$en[12]\t$en[13]\t$en[14]\t$en[15]\t$en[16]\t$en[17]\t$en[18]\t$en[19]\t$en[20]\t$en[21]\n";
#        print "$e[10]\n";
#	}
#	else {
#		print "$array[$r]\n";
#	}
	$r++;$k++;
	}

#If unequal
	else {
		print "$array[$r]\n";
	        $val=$en[10];
		$r++;$k++;
	}

}
# Last line


}



