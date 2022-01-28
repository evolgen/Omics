#!/usr/bin/bash

set -e -o pipefail

module load samtools trimal gcc ;

printf "" >testcox_mitos.cds.fa ;

count=0 ;    
cat ~/RGP/scripts/selection/names_sebastes.77.txt | sort -u |
    while read speciesname; do
        echo "  $speciesname" ;
        printf ">${speciesname}\n" >>testcox_mitos.cds.fa ;
        find /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/sequence_files/ -type f -name "cox*.FAS.aln.fa" |
            while read trimfasta; do
                count=$((count+1));
                printf "$count ";
                samtools faidx $trimfasta ${speciesname} | 
                    grep -v '^>' | tr "\n" " " | sed -e 's/ //g' >>testcox_mitos.cds.fa
            done  
        printf "\n" >>testcox_mitos.cds.fa ;
        echo ;
    done

echo "Combined all sequences"        
                    
bash ~/RGP/scripts/selection/hyphy_dnds.run.sh Vertebrate-mtDNA testcox_mitos.cds.fa /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick dnds_trees/testcox.mitos


