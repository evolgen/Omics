#!/usr/bin/bash

set -e -o pipefail

###  \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/crossprots_*.gff | grep -v Gene | parallel -j 1 sh ~/RGP/scripts/parallel_extract_introns_gff.sh

file1=$@

module load bedtools

grep -w mRNA $file1 | awk '{if($5>$4){print $1"\t"$4"\t"$5"\t"$7} else print $1"\t"$5"\t"$4"\t"$7}' | sort -k1,1 -k2,2n -k3,3nr >${file1/%.gff/.mRNA};
grep -w cds $file1 | awk '{if($5>$4){print $1"\t"$4"\t"$5} else print $1"\t"$5"\t"$4}' | sort -k1,1 -k2,2n -k3,3nr  >${file1/%.gff/.cds};
bedtools subtract -a ${file1/%.gff/.mRNA} -b ${file1/%.gff/.cds} |
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1,"ProSplign","intron",$2,$3,"1",$4,".\tmult="NR";pri=4;src=P"}' | 
    sort -k1,1 -k4,4n -k5,5n \
        >${file1/%.gff/.hints};
echo ${file1/%.gff/.hints};

