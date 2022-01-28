#!/usr/bin/bash 

set -e -o pipefail

module load samtools bedtools

mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB/acti_perspecies/;

cat /global/home/users/rohitkolora/RGP/scripts/nucleartree/list_nucleargenomes_SEB.txt |
    while read species; do
        echo "Working on ${species}";
        printf ">${species}\n" >/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB/acti_perspecies/${species}.CDS.pranktrim.FAS;
        \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB/acti_pergene/* | 
            grep -v -e reduced -e fas -e FAS -e RAxML -e fai | 
            while read acti_gene; do
                printf " ${acti_gene}" | sed -e 's/\/.*\///' ;
                samtools faidx ${acti_gene} ${species} |
                    grep -v '^>' |
                    tr "\n" " " |
                    sed -e 's/ //g' \
                        >>/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB/acti_perspecies/${species}.CDS.pranktrim.FAS;
            done
            echo "" >>/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB/acti_perspecies/${species}.CDS.pranktrim.FAS;
            printf "\nDone concatenating ${species}\n";
    done        

cat /global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB/acti_perspecies/*.FAS \
    >/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB/acti_seb.pranktrim.FASTA

