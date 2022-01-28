#!/usr/bin/bash

set -e -o pipefail

module load samtools trimal gcc ;

printf "" >combined_mitos.cds.fa ;

count=0 ;    
cat ~/RGP/scripts/selection/names_sebastes.77.txt | sort -u |
    while read speciesname; do
        echo "  $speciesname" ;
        printf ">${speciesname}\n" >>combined_mitos.cds.fa ;
        find /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/sequence_files/ -type f -name "*.FAS.aln.fa" |
            while read trimfasta; do
                count=$((count+1));
                printf "$count ";
                samtools faidx $trimfasta ;
                samtools faidx $trimfasta ${speciesname} | 
                    grep -v '^>' | tr "\n" " " | sed -e 's/ //g' >>combined_mitos.cds.fa
            done  
        printf "\n" >>combined_mitos.cds.fa ;
        echo ;
    done

echo "Combined all sequences"        

trimal -gt 1 -in combined_mitos.cds.fa -out combined_mitos.cds.nex -nexus
bash ~/RGP/scripts/selection/hyphy_dnds.run.sh Vertebrate-mtDNA combined_mitos.cds.nex /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick dnds_trees/combined.mitos



