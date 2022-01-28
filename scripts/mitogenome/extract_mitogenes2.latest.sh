#!/usr/bin/bash

set -e

module load seqtk samtools bedtools

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/MITOS2/result.bed | 
    grep -v -e melenopus -e S-ruberrimus_SEB-74 |
    head -n 1 |
    xargs cat |
    grep -v -e O[HL] -e rrn -e trn |
    sed -e '/^$/d' | 
    awk -F'\t' '{ print $4}' \
        >~/RGP/scripts/mitogenome/list_mitogenes2.txt

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.fasta | grep -v 'work' |
    grep -v -e melenopus -e S-ruberrimus_SEB-74 |
    sed -e 's/\/mitogenome\/.*//' | 
    sed -e 's/.*\/output\///' | 
    sed -e 's/\//\t/' |
    sort -k1,1 -k2,2 | 
    awk -F'\t' '{print $1"-"NR"\t"$2}' \
        >~/RGP/scripts/mitogenome/list_names2.txt

awk -F'\t' '{print $1}' ~/RGP/scripts/mitogenome/list_names2.txt \
    >~/RGP/scripts/mitogenome/list_species2.txt
awk -F'\t' '{print $2}' ~/RGP/scripts/mitogenome/list_names2.txt \
    >~/RGP/scripts/mitogenome/list_samples2.txt    

if [ -e "~/RGP/scripts/mitogenome/list_problematic2.txt" ]; then 
    rm ~/RGP/scripts/mitogenome/list_problematic2.txt; 
fi

count=1
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.fasta | 
    grep -v -e melenopus -e S-ruberrimus_SEB-74 |
    grep -v 'work' |
    fgrep -f ~/RGP/scripts/mitogenome/list_samples2.txt | 
        while read fasta1; do
            working=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*/\/mitogenome\//');
            speciesname=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' | sed -e 's/.*\/output\///' | sed -e 's/\/.*//');
            samplename=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' | sed -e 's/.*\/output\///' | sed -e 's/.*\///');

            samtools faidx $fasta1;
            size=$(awk '{print $2}' ${fasta1}.fai);

            printf "\tProcessing - ${speciesname}\n";
            mkdir -p ${working}/MITOS2/CDS;

            cat ~/RGP/scripts/mitogenome/list_mitogenes2.txt | while read genename; do
                printf "\tExtracting - ${genename}\n";
                grep -w "${genename}" ${working}/MITOS2/result.bed \
                    >${working}/MITOS2/CDS/${genename}.bed.edit

                check=$(wc -l ${working}/MITOS2/CDS/${genename}.bed.edit)
                if [ "${check}" == 0 ]; then
                    printf "\n\n\t\tERROR - ${genename} missing in ${working}\n\n"
                    exit 1
                fi

                awk -v size="$size" '{if($3<$2) {print $1"\t"$2"\t"size"\t"$1"\t"0"\t"$3} else print $1"\t"$2"\t"$3}' \
                    ${working}/MITOS2/CDS/${genename}.bed.edit >${working}/MITOS2/CDS/${genename}.bed        ### Circularity Problem
                rm ${working}/MITOS2/CDS/${genename}.bed.edit

                printf ">${speciesname} ${genename}\n" \
                    >${working}/MITOS2/CDS/${genename}.FAS
                strand=$(awk '{print $NF}' ${working}/MITOS2/CDS/${genename}.bed);
                if [ "${strand}" == "-" ]; then                 ### STRAND for REV COMPLEMENT
                    bedtools getfasta -fi ${fasta1} -bed ${working}/MITOS2/CDS/${genename}.bed |
                        grep -v '^>' |
                        rev |
                        tr "[ATGCatgc]" "[TACGtacg]" \
                            >>${working}/MITOS2/CDS/${genename}.FAS
                else 
                    bedtools getfasta -fi ${fasta1} -bed ${working}/MITOS2/CDS/${genename}.bed |
                        grep -v '^>' \
                            >>${working}/MITOS2/CDS/${genename}.FAS 
                fi    

                problem_count=$(ls -l ${working}/MITOS2/CDS/[a-z]*.FAS | awk '$5<50' | wc -l)    
                if [ "${problem_count}" -gt 0 ]; then
                    printf "${speciesname}\t${fasta1}" >>~/RGP/scripts/mitogenome/list_problematic2.txt
                fi
                count=$((count+1))
                done
        done

        if [ -e "~/RGP/scripts/mitogenome/list_problematic.txt2" ]; then 
            problems=$(wc -l ~/RGP/scripts/mitogenome/list_problematic2.txt);
            print "\n\nPLEASE CHECK FOR ${problems} ERROR SAMPLES :  - ~/RGP/scripts/mitogenome/list_problematic2.txt\n";
        fi



