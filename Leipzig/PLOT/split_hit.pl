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

if($e[0] eq $en[0] && $e[5] eq $en[5] && $e[4] eq $en[4] && $en[1] > $e[2]+25) {
#        print OUT "$array[$k]\n$array[$r]\n";
        print "$array[$k]\n$array[$r]\n";
        $r++;$k++;
}

elsif($en[1]>=$e[1] &&  $en[2]<=$e[2] && $en[6]<=$e[6]) {
        $r++;
}

else {
$k++;$r++;
}

}


=cut
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

