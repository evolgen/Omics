#!/usr/bin/bash

set -e -o pipefail

if [ ! -e "/global/home/users/rohitkolora/RGP/scripts/orthology/files_finale.list.masurca.txt" ]; then
    echo "    Create input file files_finale.list.txt"
    exit 1
fi

mkdir -p input_data;
cat /global/home/users/rohitkolora/RGP/scripts/orthology/files_finale.list.masurca.txt | grep -v FREEZE | while read file1; do 
    file2=$(echo $file1 | sed -e 's/.gff3$/.proteins.faa/'); 
    species=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/.*\/WTDBG\///' -e 's/.*\/FALCON\///' -e 's/\/.*//'); 

    printf "Creation of \t ${species}\n";
    echo $file1;
    awk -F'\t' -v species="${species}" 'BEGIN{OFS="\t"} {gsub(".cds;.*",";",$9); gsub("ID=","ID="species".",$9); print $0}' $file1 | fgrep -w CDS >input_data/${species}.gff; 
    sed -e "/^>/ s/^>/>${species}./" $file2 >input_data/${species}.faa; 
done


