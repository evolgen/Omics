#!/usr/bin/bash

set -e -o pipefail

proteinorthotsv="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/proteinortho_finale/SCOs_76/76.rockfish.poff.tsv.subset"
inputdirec="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/proteinortho_finale/input_data"
outputdir="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/proteinortho_finale/SCOs_76/fasta"

module load samtools;

count=0; 
tail -n +2 ${proteinorthotsv} | 
    cut -f4- | 
    while read orthogroup; do 
        count=$((count+1));
        printf "Group\t"${count}"\n";
        printf "" >${outputdir}/${count}.faa; 
        printf "" >${outputdir}/${count}.fna; 
        echo ${orthogroup} | tr " " "\n" | sed '/^$/d' |
        while read seqname; do 
            species=$(echo ${seqname} | sed -e 's/\..*//'); 
            printf ">${species}\n" >>${outputdir}/${count}.faa;
            printf ">${species}\n" >>${outputdir}/${count}.fna;
            samtools faidx ${inputdirec}/${species}.faa ${seqname} |
                grep -v '^>' >>${outputdir}/${count}.faa; 
            samtools faidx ${inputdirec}/${species}.fna ${seqname} |
                grep -v '^>' >>${outputdir}/${count}.fna; 
            printf "  "${species}"="${seqname}"\n"; 
        done; 

    done
printf "\nDone Extracting\n\n";


