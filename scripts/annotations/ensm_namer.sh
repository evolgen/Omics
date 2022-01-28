#!/usr/bin/bash

set -e 

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/Funannotate/predict_results/ensembl99_simhits_BRK_*.proteins.txt | 
    grep -v ruberrimus | 
    while read ensmfile; do
        dir1=$(dirname "$ensmfile"); 
        species=$(dirname "$ensmfile" | sed -e 's/\/Funannotate.*//' -e 's/.*\///') ; 
        pepfile=$(dirname "$ensmfile" | sed -e "s/$/\/Final_filt_BRK_${species}.proteins.fa.out.pepnames/")
        awk -F'\t' '$2!="*"' $ensmfile | 
            fgrep -w -f <(cut -f1 $pepfile) $ensmfile | 
            cut -f1,2 | awk -F'\t' '$2!="*"' |
            sed -e "s/^/${species}./" -e 's/\t\*\t/\tNA_UNKNOWN\t/' |
            awk -F'\t' '!x[$1]++' >${dir1}/Final_filt_ensembl99_simhits_BRK_${species}.proteins.txt.top; 
        echo ${dir1}/Final_filt_ensembl99_simhits_BRK_${species}.proteins.txt.top ;
        done 
         
#\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/Funannotate/predict_results/Final_filt_ensembl99_simhits_BRK_*.proteins.txt.top |
#    while read file1; do
#        cat $file1 |
#            while read line; do 
#                seqname=$(echo $line | cut -d' ' -f1) 
#                ensname=$(echo $line | cut -d' ' -f2); 
#                speciesname=$(echo $line | cut -d' ' -f3);
#                printf $speciesname"."$seqname"\t"; 
#                fgrep -w -m 1 -e "$ensname" /global/scratch2/rohitkolora/databases/ensembl99/vertebrates.pep.names.fish; 
#            done 
#    done >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/refspecies_ensnames.txt ;   

echo "  $species" ;

