#!/usr/bin/bash

set -e -o pipefail

module load samtools trimal gcc ;

printf "" >combined_nuclear_busco.cds.fa ;

count=0 ;    
cat ~/RGP/scripts/selection/names_sebastes.77.txt | sort -u |
    while read speciesname; do
        echo "  $speciesname" ;
        printf ">${speciesname}\n" >>combined_nuclear_busco.cds.fa ;
        find /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/nucleargenome/BUSCO/sequence_files/*/ -type f -name 'NT.FNA.trim.fas' |
            while read trimfasta; do
                count=$((count+1));
                printf "$count ";
                samtools faidx $trimfasta ${speciesname} | 
                    grep -v '^>' | tr "\n" " " | sed -e 's/ //g' >>combined_nuclear_busco.cds.fa ;
            done  
        printf "\n" >>combined_nuclear_busco.cds.fa ;
        echo ;
    done

echo "Combined all sequences"        
                    

