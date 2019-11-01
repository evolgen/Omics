#!/usr/bin/bash

set -e

module load seqtk samtools bedtools


\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.fasta | 
    grep -v 'work' |
    fgrep -f ~/RGP/scripts/mitogenome/list_samples2.txt | 
        while read fasta1; do
            working=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*/\/mitogenome\//');
            speciesname=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' | sed -e 's/.*\/output\///' | sed -e 's/\/.*//');
            samplename=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' | sed -e 's/.*\/output\///' | sed -e 's/.*\///');

            touch "${working}/MITOS2/CDS/ALL.FAS"; ####    if [ ! -e "${working}/MITOS2/CDS/ALL.FAS" ]; then
            samtools faidx $fasta1;

            printf "\tProcessing - ${speciesname}\n";

            printf ">${speciesname}.${samplename}\n" \
                        >${working}/MITOS2/CDS/ALL.FAS

            cat ~/RGP/scripts/mitogenome/list_mitogenes2.txt | while read genename; do
                printf "\t${genename} ";
                cat ${working}/MITOS2/CDS/${genename}.FAS |
                    grep -v '^>' |
                    sed -e '/^$/d' |
                    tr "\n" " " |
                    sed -e 's/ //g' \
                        >>${working}/MITOS2/CDS/ALL.FAS

            done
            echo "";

        done
ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/MITOS2/CDS/ALL.FAS | 
    xargs cat | 
    sed -e 's/>/\n>/' | 
    sed -e '/^$/d' \
        >CDS.all.FAS

