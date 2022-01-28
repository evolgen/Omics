#!/usr/bin/bash

set -e -o pipefail

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/crossprots_*.gff | \
    grep -v Genes | \
    while read line; do \
        if [ ! -e "${line}3" ]; then \
            echo $(basename "$line");
            genome=$(dirname "$line" | sed -e 's/$/\/final.genome.scf.FAS/'); 
            if [ ! -e "${genome}".fai ]; then
                samtools faidx ${genome};    
            fi    
                grep -v '^#' ${line} | \
                perl ~/RGP/scripts/annotations/fix_genomecoords_gfforgtf.pl ${genome}.fai - \
                    >${line}.fix;
            echo "${line}.fix";
            /global/scratch2/rohitkolora/Software/gffread/gffread ${line}.fix \
                -g ${genome} --adj-stop \
                -V -C -M -K -Z -F \
                -o ${line}3 1>${line/%.gff/.err} 2>&1; 
        fi; 
    done


