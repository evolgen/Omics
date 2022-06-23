#!/usr/bin/perl

# perl blast_threshold.pl infile > outfile

$filename = $ARGV[0];
open($filehandle, '<', $filename) or die "Could not open $filename\n";
@arr=<$filehandle>;
chomp(@arr);
close(F);

local $" = "\t";			#for printing tab in array

for ($j=0;$j<=$#arr;$j++)
{
@blast=split(/\t/,$arr[$j]);
if ($blast[10] <= 0.000001)
{
print "@blast\n";
}
}



