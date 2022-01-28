#!/usr/bin/bash

set -e

printf "" > /global/scratch2/rohitkolora/Rockfish/Genomes/rearrangements/rockfish.per_qry.by_aleu.names ;
printf "" > /global/scratch2/rohitkolora/Rockfish/Genomes/rearrangements/rockfish.per_qry.by_aleu.RC ;

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/Pairwise/SYRI/*/syri1.Sebastes_aleutianus.*.chrnames | 
    grep -v Medaka | 
    while read file1; do 
        file2=$(dirname "$file1" | sed -e 's/$/\/log_syri1_rcnames/'); 
        species=$(basename "$file1" | sed -e 's/\.chrnames//' -e 's/.*\.//'); 
        echo $species ;
        grep -A 100 -e '###QRY-names###' $file1 | grep -v -e PGA -e '#' | sed -e "s/^/${species}\t/" -e 's/>//g' | 
            awk -F'\t' '$3!="" {print $1"\t"$3"\t"$2}' >>/global/scratch2/rohitkolora/Rockfish/Genomes/rearrangements/rockfish.per_qry.by_aleu.names ; 
        sed -e "s/^/${species}\t/" -e 's/$/\tRC/' $file2 | sort -k2,2V >>/global/scratch2/rohitkolora/Rockfish/Genomes/rearrangements/rockfish.per_qry.by_aleu.RC ; 
    done

