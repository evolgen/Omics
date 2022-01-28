#!/usr/bin/bash

samfile=$1
bamfile=$2

if [ "$#" -ne 2 ]; then
    echo "   Need two files - sam and bam";
    exit 1;
fi    

module load samtools

awk -F'\t' 'BEGIN{OFS="\t"} {gsub("M","=",$6); print $0}' $samfile |
    samtools sort -@ 20 -o $bamfile - ;
    

