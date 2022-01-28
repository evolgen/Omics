#!/usr/bin/bash

set -e

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/INTERPROSCAN_Fun/interproscan.tsv | while read file1; do 
    species=$(echo $file1 | sed -e 's/.*FREEZE\///' -e 's/\/INTERPRO.*//' -e 's/.*\///'); 
    file2=$(echo $file1 | sed -e "s/INTERPROSCAN_Fun\/interproscan.tsv$/Funannotate\/predict_results\/BRK_${species}.gff3/"); 
    cat ~/RGP/scripts/annotations/visualopsin.panther.txt | while read opsindet; do 
        name=$(echo $opsindet | sed -e 's/.* //' -e 's/.*\t//'); 
        opsin=$(echo $opsindet | sed -e 's/ .*//' -e 's/\t.*//' ); 
        printf ${species}"\t"$name"\t"; 
        fgrep -w "${opsin}" $file1 | cut -f1,5 | sed -e 's/\..*\t/\t/' | sort -k2,2V | uniq -D -f1 | sort -u | cut -f1 | while read line; do 
            grep -m 1 -w -e "$line" $file2; 
        done | 
        bedtools cluster -i - | wc -l; 
    done; 
done 

