#!/usr/bin/perl
$file = $ARGV[0];
open(Input,"<","$file");

@array=<Input>;
$k=0;

$ary_len=$#array;
while($array[$k])
{
chomp $array[$k];
chomp $array[$k+1];
@e=split /\t/,$array[$k];
@en=split /\t/,$array[$k+1];
if($k==$ary_len && $e[2]!=$e[3])
{
$gap=$e[3]-$e[2];
print "$e[0]\t$e[2]\t$e[3]\t$gap\t$e[4]\n"; }

if($k==0 && $e[1]!=0)
{ print "$e[0]\t0\t$e[1]\t$e[1]\t$e[4]\n"; }

if($e[0] eq $en[0] )
        {
        $gap=$en[1]-$e[2];
                print "$e[0]\t$e[2]\t$en[1]\t$gap\t$e[4]\n";
        }
$k++;
}

#126755	bil_Backbone_1969	504803	1520	135680	-	vir_Backbone_2	146312	0	136028
#0	1			2	3	4	5	6		7	8	9
#0	9			10	11	12	8	13		14	15	16
#-	0			3	1	2	5	-		-	-	-
#-	0			3	1	2	4	-		-	-	-
#bil_Backbone_10	0	894	894	

