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

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/SEB"
mkdir -p ${working}/${group}_pergene

#if [ ! -e "~/RGP/scripts/nucleartree/list_nucleargenomes_SEB.txt" ]; then
#printf "Creating genomes list \n";
#grep -e Sebastes_ -e Sebastolobus_alascanus ~/RGP/scripts/nucleartree/list_nucleargenomes_acti.txt |
#    grep -v -e Sebastes_vulpes -e Sebastes_hopkinsi -e Sebastes_borealis -e Sebastes_cheni -e Sebastes_melanopus \
#        >~/RGP/scripts/nucleartree/list_nucleargenomes_SEB.txt
#fi        

number=$(wc -l ~/RGP/scripts/nucleartree/list_nucleargenomes_SEB.txt | awk '{print $1}');
printf " Looking for atleast ${number} Species\n";

#if [ ! -e "~/RGP/scripts/nucleartree/list_scp_${group}_SEB-files.txt" ]; then
#printf "\tCreating file list\n";
#\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/*.fna |
#        fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes_SEB.txt |
#                sed -e 's/\/.*\///' >~/RGP/scripts/nucleartree/list_scp_${group}_SEB-files.txt;
#fi                

#if [ ! -e "~/RGP/scripts/nucleartree/list_scp_${group}_SEB-genes.txt" ];  then
#printf "\tCreating gene list\n";
#        cat ~/RGP/scripts/nucleartree/list_scp_${group}_SEB-files.txt | 
#        sort | uniq -c |
#        awk -v number="$number" '$1==number' |
#        sed -e '/^$/d' \
#            >~/RGP/scripts/nucleartree/list_scp_${group}_SEB-genes.txt
#fi            

    awk '{print $2}' ~/RGP/scripts/nucleartree/list_scp_${group}_SEB-genes.txt \
        >~/RGP/scripts/nucleartree/list_scp_${group}_SEB-genes.filt.txt


count=1
cat ~/RGP/scripts/nucleartree/list_scp_${group}_SEB-genes.filt.txt | 
    while read gene; do
        fastaname=$(echo ${gene} | sed -e 's/.*\///' -e 's/\.fna//');    
        printf "Executing for ${fastaname} \n";
        touch ${working}/${group}_pergene/${fastaname}.FAS;

        \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/ |
            fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes_SEB.txt |
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
    parallel -j 20 sh ~/RGP/scripts/nucleartree/run_raxml_pergene_SEB.sh 


