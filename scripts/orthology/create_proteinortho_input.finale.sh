#!/usr/bin/bash

set -e -o pipefail

if [ ! -e "files_finale.list.txt" ]; then
    echo "    Create input file files_finale.list.txt"
    exit 1
fi

mkdir -p input_data;
cat files_finale.list.txt | grep -v FREEZE | while read file1; do 
    file2=$(echo $file1 | sed -e 's/finale_seq.CDS.gff/finale_seq.aa/'); 
    species=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/.*\/WTDBG\///' -e 's/.*\/FALCON\///' -e 's/\/.*//'); 

    printf "Creation of \t ${species}\n";
    echo $file1;
#    awk -F'\t' -v species="${species}" 'BEGIN{OFS="\t"} {gsub("\";.*",";",$0); gsub("transcript_id.*\"g","ID="species".g",$0); print $0}' $file1 >input_data/${species}.gff; 
    cp $file1 input_data/${species}.gff;
    cp $file2 input_data/${species}.faa; 
done

cat files_finale.list.txt | grep -e FREEZE | while read file1; do
    file2=$(echo $file1 | sed -e 's/crossprot2aug.norepeats.CDS.gff/Annot_crossprot.faa/');
    species=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/.*\/WTDBG\///' -e 's/.*\/FALCON\///' -e 's/\/.*//');
    printf "Creation of \t ${species}\n";
    echo $file1;
       awk -F'\t' -v species="${species}" 'BEGIN{OFS="\t"} {gsub("\";.*",";",$0); gsub("transcript_id.*\"g","ID="species".g",$0); print $0}' $file1 >input_data/${species}.gff;
       sed -e '/^>/ s/>FALCON_/>/' -e 's/>WTDBG_/>/' $file2 >input_data/${species}.faa;
done

sed -e 's/Sebastes_ruberrimus/Sebastes_miniatus/g' /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_ruberrimus/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gff | awk -F'\t' -v species="Sebastes_miniatus" 'BEGIN{OFS="\t"} {gsub("\";.*",";",$0); gsub("transcript_id.*\"g","ID="species".g",$0); print $0}' >input_data/Sebastes_miniatus.gff
sed -e 's/Sebastes_ruberrimus/Sebastes_miniatus/g' -e '/^>/ s/>FALCON_/>/' /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_ruberrimus/FILTERED/CROSSPROT/Annot_crossprot.faa >input_data/Sebastes_miniatus.faa


