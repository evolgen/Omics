#!/usr/bin/bash

set -e

module load seqtk samtools

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/MITOS/result.bed | 
    head -n 1 |
    xargs cat |
    grep -v -e O[HL] -e rrn -e trn |
    sed -e '/^$/d' | 
    awk -F'\t' '{ print $4}' >~/RGP/scripts/mitogenome/list_mitogenes.txt

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.result/*.fasta | grep -v 'work' |
    sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/\//\t/' |
    sort -k1,1 -k2,2 | awk -F'\t' '!x[$1]++ {print $1"-"NR"\t"$2}' >~/RGP/scripts/mitogenome/list_names.txt

awk -F'\t' '{print $1}' ~/RGP/scripts/mitogenome/list_names.txt >~/RGP/scripts/mitogenome/list_species.txt
awk -F'\t' '{print $2}' ~/RGP/scripts/mitogenome/list_names.txt >~/RGP/scripts/mitogenome/list_samples.txt    

if [ -e "~/RGP/scripts/mitogenome/list_problematic.txt" ]; then 
    rm ~/RGP/scripts/mitogenome/list_problematic.txt; 
fi

count=1
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.result/*.fasta | 
    grep -v 'work' |
#    fgrep -f ~/RGP/scripts/mitogenome/list_samples.txt | 
        while read fasta1; do
            working=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*/\/mitogenome\//');
            speciesname=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/\/.*//');
            samplename=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/.*\///');

            if [ ! -e "${working}/MITOS/CDS/ALL.FAS" ]; then
            samtools faidx $fasta1;
            size=$(awk '{print $2}' ${fasta1}.fai);

            printf "\tProcessing - ${speciesname}\n";
            mkdir -p ${working}/MITOS/CDS;
            touch ${working}/MITOS/CDS/${speciesname}.fas
            cat ~/RGP/scripts/mitogenome/list_mitogenes.txt | while read genename; do
                printf "\tExtracting - ${genename}\n";
               grep -w "${genename}" ${working}/MITOS/result.bed >${working}/MITOS/CDS/${genename}.bed.edit
                awk -v size="$size" '{if($3<$2) {print $1"\t"$2"\t"size"\t"$1"\t"0"\t"$3} else print $1"\t"$2"\t"$3}' \
                    ${working}/MITOS/CDS/${genename}.bed.edit >${working}/MITOS/CDS/${genename}.bed        ### Circularity Problem
                rm ${working}/MITOS/CDS/${genename}.bed.edit
                seqtk subseq ${fasta1} ${working}/MITOS/CDS/${genename}.bed |
                    grep -v '^>' |
                    sed -e "1i >${speciesname} ${genename}" >${working}/MITOS/CDS/${genename}.FAS
                cat ${working}/MITOS/CDS/${genename}.FAS | grep -v '^>' >>${working}/MITOS/CDS/${speciesname}.fas    
            done
            printf ">${speciesname}-${count}\n" >${working}/MITOS/CDS/ALL.FAS
###            printf ">${speciesname}.${samplename}\n" >${working}/MITOS/CDS/ALL.FASTA
###            cat ${working}/MITOS/CDS/ALL.FAS | grep -v '^>' | tr "\n" " " | sed -e 's/ //g' -e 's/>/\n>/' -e '/^$/d'  >>${working}/MITOS/CDS/ALL.FASTA
            cat ${working}/MITOS/CDS/${speciesname}.fas | 
                tr "\n" " " | 
                sed -e 's/ //g' -e 's/>/\n>/' -e '/^$/d'  >>${working}/MITOS/CDS/ALL.FAS
                rm ${working}/MITOS/CDS/${speciesname}.fas
            fi    
            problem_count=$(ls -l ${working}/MITOS/CDS/[a-z]*.FAS | awk '$5<50' | wc -l)    
            if [ "${problem_count}" -gt 0 ]; then
                printf "${speciesname}\t${fasta1}" >>~/RGP/scripts/mitogenome/list_problematic.txt
            fi
            count=$((count+1))
        done

        if [ -e "~/RGP/scripts/mitogenome/list_problematic.txt" ]; then 
            problems=$(wc -l ~/RGP/scripts/mitogenome/list_problematic.txt);
            print "\n\nPLEASE CHECK FOR ${problems} ERROR SAMPLES :  - ~/RGP/scripts/mitogenome/list_problematic.txt\n";
        fi



