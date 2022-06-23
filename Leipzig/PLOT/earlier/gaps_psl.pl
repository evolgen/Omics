#!/usr/bin/perl
$file = $ARGV[0];
open(Input,"<","$file");
#open (OUT,">","$file.gaps");
@array=<Input>;
$k=0;
#print OUT "#query\tstart_gap\tend_gap\tleng_gap\tsign\n";
$ary_len=$#array;
while($array[$k])
{
chomp $array[$k];
chomp $array[$k+1];
@e=split /\t/,$array[$k];
@en=split /\t/,$array[$k+1];
if($k==$ary_len && $e[12]!=$e[10])
{
$gap=$e[10]-$e[12];
print "$e[9]\t$e[12]\t$e[10]\t$gap\t$e[8]\n"; }

if($k==0 && $e[11]!=0)
{ print "$e[9]\t0\t$e[11]\t$e[11]\t$e[8]\n"; }

if($e[1] eq $en[1] )
        {
        $gap=$en[11]-$e[12];
                print "$e[9]\t$e[12]\t$en[11]\t$gap\t$e[8]\n";
        }
$k++;
}

#126755	bil_Backbone_1969	504803	1520	135680	-	vir_Backbone_2	146312	0	136028
#0	1			2	3	4	5	6		7	8	9
#0	9			10	11	12	8	13		14	15	16


