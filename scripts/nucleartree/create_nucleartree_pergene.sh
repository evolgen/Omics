#!/usr/bin/bash

set -e -o pipefail

group=$1

if [ "$#" -ne 1 ]; then
    printf "Illegal number of parameters:\t     $#      \n";
    exit 1;
fi        

if [ "$group" != "acti" ] && [ "$group" != "vert" ] && [ "$group" != "cvg" ]; then
    printf "Illegal BUSCO-SCG database - need either acti or vert or cvg\n\n";
fi    

if [ "$group" == "verti" ]; then
    limit=2000
fi    

if [ "$group" == "acti" ]; then
    limit=3400
fi

if [ "$group" == "cvg" ]; then
    limit=150
fi

module load seqtk samtools bedtools

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree"
mkdir -p ${working}/${group}_pergene

#if [ ! -e "~/RGP/scripts/nucleartree/list_nucleargenomes_${group}.txt" ]; then
#    printf "\tDoes not exit : ~/RGP/scripts/nucleartree/list_nucleargenomes_${group}.txt\n\n";
#    exit 1
#fi
#if [ ! -e "~/RGP/scripts/nucleartree/list_scp_${group}-genes.filt.txt" ]; then
#    printf "\tDoes not exit : ~/RGP/scripts/nucleartree/list_scp_${group}-genes.filt.txt\n\n";
#    exit 1
#fi


count=1
cat ~/RGP/scripts/nucleartree/list_scp_${group}-genes.filt.txt | 
    while read gene; do
        fastaname=$(echo ${gene} | sed -e 's/.*\///' -e 's/\.fna//');    
        printf "Executing for ${fastaname} \n";
        touch ${working}/${group}_pergene/${fastaname}.FAS;

        \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/ |
            fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes_${group}.txt |
            while read folder; do
                speciesname=$(echo ${folder} | sed -e 's/.*\/output\///' -e 's/\/.*//');
                printf " Adding ${speciesname} \n\t";


                if [ ! -e "${folder}/${gene}" ]; then
                    printf "ERROR : Missing in ${speciesname}\n ${folder}/${gene}\n\n";
                    exit 1
                fi
                
                printf ">${speciesname}\n" \
                    >>${working}/${group}_pergene/${fastaname}.FAS
                cat ${folder}/${gene} |
                    grep -v '^>' |
                    sed -e 's/taa$//' -e 's/tag$//' -e 's/tga$//' \
                        >>${working}/${group}_pergene/${fastaname}.FAS 

                printf "\n" >>${working}/${group}_pergene/${fastaname}.FAS
                done

                printf " Finished with ${fastaname}\n";
            done


\ls -d ${working}/${group}_pergene/*.FAS |
    parallel -j 32 sh ~/RGP/scripts/nucleartree/run_raxml_pergene.sh 


