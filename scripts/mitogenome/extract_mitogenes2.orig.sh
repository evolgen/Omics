#!/usr/bin/bash

set -e

module load seqtk samtools bedtools

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/MITOS2/result.bed | 
    head -n 1 |
    xargs cat |
    grep -v -e O[HL] -e rrn -e trn |
    sed -e '/^$/d' | 
    awk -F'\t' '{ print $4}' \
        >~/RGP/scripts/mitogenome/list_mitogenes2.txt

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.fasta | grep -v 'work' |
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
    grep -v 'work' |
    fgrep -f ~/RGP/scripts/mitogenome/list_samples2.txt | 
    grep -e Ade | 
        while read fasta1; do
            working=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*/\/mitogenome\//');
            speciesname=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' | sed -e 's/.*\/output\///' | sed -e 's/\/.*//');
            samplename=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' | sed -e 's/.*\/output\///' | sed -e 's/.*\///');

            touch "${working}/MITOS2/CDS/ALL.FAS"; ####    if [ ! -e "${working}/MITOS2/CDS/ALL.FAS" ]; then
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
##                awk 'OFS="\t" {if($2>$3) print $1, $3, $2, ".", ".", "-"; else print $0, ".", ".", "+"}'  ${working}/MITOS2/CDS/${genename}.bed | 

                printf ">${speciesname} ${genename}\n" \
                    >${working}/MITOS2/CDS/${genename}.FAS
                strand=$(awk '{print $NF}' ${working}/MITOS2/CDS/${genename}.bed);
                if [ "${strand}" == "-" ]; then                 ### STRAND for REV COMPLEMENT
                    bedtools getfasta -fi ${fasta1} -bed ${working}/MITOS2/CDS/${genename}.bed |
##                seqtk subseq ${fasta1} ${working}/MITOS2/CDS/${genename}.bed |
                        grep -v '^>' |
                        rev |
                        tr "[ATGCatgc]" "[TACGtacg]" \
                            >>${working}/MITOS2/CDS/${genename}.FAS
                else 
                    bedtools getfasta -fi ${fasta1} -bed ${working}/MITOS2/CDS/${genename}.bed |
                        grep -v '^>' \
                            >>${working}/MITOS2/CDS/${genename}.FAS 
                fi    
####                sed -e "1i >${speciesname} ${genename}" >${working}/MITOS2/CDS/${genename}.FAS

                    cat ${working}/MITOS2/CDS/${genename}.FAS | 
                        grep -v '^>' \
                            >>${working}/MITOS2/CDS/${speciesname}.fas    

                    printf ">${speciesname}.${samplename}\n" \
                        >${working}/MITOS2/CDS/ALL.FAS
                    cat ${working}/MITOS2/CDS/${genename}.FAS |
                        grep -v '^>' |
                        sed -e '/^$/d' \
                            >>${working}/MITOS2/CDS/ALL.FAS

                    ####    else            
####        printf "Exists  -   ${working}/MITOS2/CDS/ALL.FAS\n";
####    fi        
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



