#!/usr/bin/bash

set -e -o pipefail

file1=$1
fileoutsum=$2
fileouttop=$3

if [ "$#" -ne 3 ]; then
    echo "   Need an input paf file and output file - sum and top";
    exit 1;
fi

columncheck=$(awk -F'\t' '{print NF}' $file1 | sort -u);
columncount=$(awk -F'\t' '{print NF}' $file1 | sort -u | wc -l);

if [ "$columncheck" -eq 16 ] ; then
    cut -f1-12 $file1 >${file1}.tmp ;
    mv ${file1}.tmp $file1 ;
    columncheck=$(awk -F'\t' '{print NF}' $file1 | sort -u) ;
fi

if [ "$columncheck" -ne 12 ] || [ "$columncount" -ne 1 ] ; then
    echo "   Strictly need a 12-column paf file validity";
    echo "   You have $columncount column types and column numbers = $columncheck";
    exit 1;
fi

printf "\t Running a valid paf file with 12 columns - $file1\n";

awk -F'\t' '$12==60 { pairsum[$1"\t"$6] += $10; pairlength[$1"\t"$6] += $11 } \
    END {for (pair in pairsum) \
    print pair"\t"pairsum[pair]"\t"pairlength[pair]"\t"pairsum[pair]/pairlength[pair] }' \
    $file1 >$fileoutsum;

sort -t$'\t' -k2,2V -k3,3nr -k4,4nr $fileoutsum | awk -F'\t' '!x[$2]++' >$fileouttop ;

