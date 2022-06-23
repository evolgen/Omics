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
if($k==$ary_len && $e[4]!=$e[2])
{
$gap=$e[2]-$e[4];
print "$e[1]\t$e[4]\t$e[2]\t$gap\t$e[5]\n"; }

if($k==0 && $e[3]!=0)
{ print "$e[1]\t0\t$e[3]\t$e[3]\t$e[5]\n"; }

if($e[1] eq $en[1] )
        {
        $gap=$en[3]-$e[4];
                print "$e[1]\t$e[4]\t$en[3]\t$gap\t$e[5]\n";
        }
$k++;
}


