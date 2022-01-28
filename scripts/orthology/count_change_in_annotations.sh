#!/usr/bin/bash

set -o -e pipefail

printf "species\tsample\tbraker\tstartstop\tcrossprot\tstep1\tafterrpm\n" >sample_annotation_count.list; 

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/ | 
    while read line; do 
        species=$(echo $line | sed -e 's/.*\/output\///' -e 's/\/.*//'); 
        sample=$(echo $line | sed -e 's/\/assembly\/.*//' -e 's/.*\/output\///' -e 's/.*\///'); 
        if [ -e "${line}/FILTERED/CROSSPROT/Annot_crossprot.faa.fai" ]; then 
            printf ${species}"\t"${sample}"\t"$(cat ${line}/BRAKER/augustus.ab_initio.aa.fai | wc -l)"\t"$(cat ${line}/FILTERED/metstop4aug.list | wc -l)"\t"$(cat ${line}/FILTERED/crossprot2aug.list | wc -l)"\t"$(cat ${line}/FILTERED/Annot_final.aa.fai | wc -l)"\t"$(cat ${line}/FILTERED/CROSSPROT/Annot_crossprot.faa.fai | wc -l)"\n"; 
        fi; 
done >>sample_annotation_count.list

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/ | while read line; do 
    species=$(echo $line | sed -e 's/.*\/FALCON\///' -e 's/.*\/WTDBG\///' -e 's/\/.*//'); 
    sample=$(echo $line | sed -e 's/.*\/FREEZE\///' -e 's/\/.*//'); 
    if [ -e "${line}/FILTERED/CROSSPROT/Annot_crossprot.faa.fai" ]; then 
        printf ${species}"\t"${sample}"\t"$(cat ${line}/BRAKER/augustus.ab_initio.aa.fai | wc -l)"\t"$(cat ${line}/FILTERED/metstop4aug.list | wc -l)"\t"$(cat ${line}/FILTERED/crossprot2aug.list | wc -l)"\t"$(cat ${line}/FILTERED/Annot_final.aa.fai | wc -l)"\t"$(cat ${line}/FILTERED/CROSSPROT/Annot_crossprot.faa.fai | wc -l)"\n"; 
    fi; 
done >>sample_annotation_count.list


