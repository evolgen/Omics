#!/usr/bin/bash

set -e -o pipefail

proteinorthotsv="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/rockfish_refgen.poff.filled.tsv.fna"
inputdirec="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/input_all_lifted"
outputdir="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Ortho_groups_fna"
dict_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/proteinnames_all_withlifted.noscaff.txt"

module load samtools;

count=0; 
tail -n +2 ${proteinorthotsv} | 
    cut -f4- | 
    while read orthogroup; do 
        count=$((count+1));
        count=$(printf "%05d\n" $count);
        printf "Group\t"${count}"\n";
        printf "" >${outputdir}/OrthoGroup.${count}.faa; 
        printf "" >${outputdir}/OrthoGroup.${count}.fna; 
        echo ${orthogroup} | tr " " "\n" | sed '/^$/d' |
        while read seqname; do 
            presence=$(fgrep -m 1 "$seqname" $dict_dir | wc -l) ;
            species=$(echo ${seqname} | sed -e 's/\..*//'); 
            if [ "$presence" -eq "1" ]; then
                samtools faidx ${inputdirec}/${species}.faa ${seqname} >>${outputdir}/OrthoGroup.${count}.faa ;
                samtools faidx ${inputdirec}/${species}.fna ${seqname} >>${outputdir}/OrthoGroup.${count}.fna ;
                printf "  "${species}"="${seqname}"\n"; 
            fi
        done; 

    done
printf "\nDone Extracting\n\n";


