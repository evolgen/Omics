#!/usr/bin/bash

set -e

cd /global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation ;
mkdir -p iterations ;

for count in `seq 1 100`; do
    mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count ;
    head -n 12 /global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/master_alignment_convervation_AA.py >/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;
    printf "\nlong_lived_list = [" >>/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;
    shuf -n 7 /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/species_list.55.txt >/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/long_lived.list ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/long_lived.list | sed -e 's/^/"/g' -e 's/$/"/' | tr "\n" "," | sed -e 's/,$//' >>/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;
    printf "]\n" >>/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;

    printf "\nshort_lived_list = [" >>/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;
    fgrep -v -f /global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/long_lived.list /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/species_list.55.txt |
        shuf -n 11 >/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/short_lived.list ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/short_lived.list |    
        sed -e 's/^/"/g' -e 's/$/"/' | tr "\n" "," | sed -e 's/,$//' >>/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;
    printf "]\n\n" >>/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;
    tail -n +23 /global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/master_alignment_convervation_AA.py >>/global/scratch2/rohitkolora/Rockfish/Genomes/RERconverge/conservation/iterations/$count/script.py ;

    printf "$count ";
done    
echo ;


