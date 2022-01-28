#!/usr/bin/bash

set -e 

printf "Species\tsample\thet_min\thet_max\n" >heterozygosit_estimates.txt; 

find /global/scratch2/rohitkolora/Rockfish/Genomes/estimates/kmer_histo -type f -name "summary.txt" | 
    while read file; do 
        name=$(dirname "$file" | sed -e 's/.*\/kmer_histo//');
         species=$(dirname "$file" | sed -e 's/.*\/kmer_histo\///' -e 's/\/.*//'); 
         sample=$(dirname "$file" | sed -e 's/.*\///'); 
         printf "$speci es\t$sample\tHeterozygosity\t"; 
         grep -w 'Heterozygosity' $file | sed -e 's/%$//g' | awk '{print $2"\t"$3}'; done >>heterozygosit_estimates.txt



