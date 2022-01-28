#!/usr/bin/bash

set -e -o pipefail

module load samtools

if [ ! -e "/global/home/users/rohitkolora/RGP/scripts/orthology/files_latest.list.masurca.txt" ]; then
    echo "    Create input file files_finale.list.txt"
    exit 1
fi

mkdir -p input_data;
cat /global/home/users/rohitkolora/RGP/scripts/orthology/files_latest.list.masurca.txt | grep -v FREEZE | while read file1; do 
    file2=$(echo $file1 | sed -e 's/\/Filt_/\//' -e 's/.proteins.faa/.gff3/'); 
    filecds=$(echo $file2 | sed -e 's/.gff3/.cds-transcripts.fa/');
    filecdsout=$(echo $file1 | sed -e 's/.proteins.faa/.cds-transcripts.fa/');
    species=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/.*\/WTDBG\///' -e 's/.*\/FALCON\///' -e 's/\/.*//'); 

    printf "Creation of \t ${species}\n";
    echo $file1;
    grep -w -e CDS $file2 | sed -e 's/\.cds;.*/;/' |
        awk -F'\t' -v species="${species}" 'BEGIN{OFS="\t"} {gsub("ID=","ID="species".",$0); print $0}' >input_data/${species}.gff; 
    sed -e "/^>/ s/>/>${species}./" $file1 >input_data/${species}.faa; 
    grep '^>' $file1 |
        sed -e 's/^>//' -e "s/${species}\.//" |
        xargs samtools faidx $filecds >$filecdsout ;
    sed -e "/^>/ s/>/>${species}./" $filecdsout >input_data/${species}.fna;    
done

