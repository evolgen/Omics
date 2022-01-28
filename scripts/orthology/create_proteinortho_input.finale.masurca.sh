#!/usr/bin/bash

set -e -o pipefail

if [ ! -e "files_finale.list.masurca.txt" ]; then
    echo "    Create input file files_finale.list.masurca.txt"
    exit 1
fi

mkdir -p input_data;
cat files_finale.list.masurca.txt | grep -v -e 'S-ruberrimus_SEB-74' | \
    while read file1; do 

    file2=$(echo $file1 | sed -e 's/finale_seq.CDS.gff/finale_seq.aa/'); 
    species=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/.*\/WTDBG\///' -e 's/.*\/FALCON\///' -e 's/\/.*//'); 

    printf "Creation of \t ${species}\n";
    echo $file1;
#    awk -F'\t' -v species="${species}" 'BEGIN{OFS="\t"} {gsub("\";.*",";",$0); gsub("transcript_id.*\"g","ID="species".g",$0); print $0}' $file1 >input_data/${species}.gff; 
    cp $file1 input_data/${species}.gff;
    cp $file2 input_data/${species}.faa; 

done

sed -e 's/Sebastes_ruberrimus/Sebastes_miniatus/g' /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/Sebastes_ruberrimus/S-ruberrimus_SEB-74/assembly/masurca/FILTERED/CROSSPROT/FINALE/finale_seq.CDS.gff | awk -F'\t' -v species="Sebastes_miniatus" 'BEGIN{OFS="\t"} {gsub("\";.*",";",$0); gsub("transcript_id.*\"g","ID="species".g",$0); print $0}' >input_data/Sebastes_miniatus.gff
sed -e 's/Sebastes_ruberrimus/Sebastes_miniatus/g' /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/Sebastes_ruberrimus/S-ruberrimus_SEB-74/assembly/masurca/FILTERED/CROSSPROT/FINALE/finale_seq.aa >input_data/Sebastes_miniatus.faa


