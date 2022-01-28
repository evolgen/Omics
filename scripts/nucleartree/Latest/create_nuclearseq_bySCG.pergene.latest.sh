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

if [ "$group" == "vert" ]; then
    limit=2000
fi    

if [ "$group" == "acti" ]; then
    limit=3400
fi

if [ "$group" == "cvg" ]; then
    limit=150
fi

module load seqtk samtools bedtools

if [ "${limit}" == "" ]; then
    printf "\n\tProblem with the minimum limit of ${group} busco SCOs - ${limit}\n\n";
    exit 0;
fi    

mkdir -p ~/RGP/scripts/nucleartree/Latest;
printf "\t Minimal limit of SCOs is ${limit}\n";

filename="~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt"
#if [ ! -f "~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt" ] ; then
#    echo "Missing - $filename"
#    exit 1
#fi

number=$(wc -l ~/RGP/scripts/nucleartree/Latest/list_nucleargenomes_${group}.uniq.latest.txt | awk '{print $1}');
printf " Creating the list of genes for ${number} unique genomes\n";
    
filename="~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt"
#if [ ! -f "${filename}" ] ; then
#    echo "Missing - $filename"
#    exit 1
#fi

scp_genes=$(cat ~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt | wc -l);

printf "\tSCO genes in total\t${scp_genes}\n";
minimal=$(($limit/10));

if [ "$scp_genes" -lt "${minimal}" ]; then
    printf "\n\tToo Few SCO genes - ${scp_genes}\n\n";
    exit 0;
fi    

printf " Starting the extraction of each gene seperately for all species\n";

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/BUSCO"
mkdir -p ${working}/${group}_per_gene;
filename="${working}/${group}_per_gene/DONE.txt"
if [ -f "${filename}" ]; then
    printf "Remove the file below to rerun extraction \n  ${filename}"
    exit 1
fi    

count=1
cat ~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt | 
    while read gene; do
        genename=$(echo ${gene} | sed -e 's/.*\///' -e 's/\.fna//');
        printf "" >${working}/${group}_per_gene/${genename}.FNA;
        printf "" >${working}/${group}_per_gene/${genename}.FAA;
        printf "${count} : ${genename}\t-\t";

        \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/ |
        fgrep -w -f ~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt |
        while read folder; do
            speciesname=$(echo ${folder} | sed -e 's/.*\/output\///' -e 's/\/.*//');

            if [ ! -e "${folder}/${genename}.fna" -o ! -e "${folder}/${genename}.faa" ]; then
                printf "ERROR : Missing in ${speciesname}\n ${folder}/${genename}.fna\n\n";
                exit 1
            fi
            
            printf " ${speciesname} ";
            printf ">${speciesname}\n" \
                >>${working}/${group}_per_gene/${genename}.FNA
            printf ">${speciesname}\n" \
                >>${working}/${group}_per_gene/${genename}.FAA
            cat ${folder}/${genename}.fna |
                grep -v '^>' | grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
                fold -w 3 | sed -e '$s/taa$//' -e '$s/tag$//' -e '$s/tga$//' |      # Stop codon craziness
                sed -e '$s/TAA$//' -e '$s/TAG$//' -e '$s/TGA$//' |
                tr "\n" " " | sed -e 's/ //g' -e '/^$/d' | 
                sed -e 's/$/\n/' >>${working}/${group}_per_gene/${genename}.FNA
            
            cat ${folder}/${genename}.faa |
                grep -v '^>' | tr "\n" " " | sed -e 's/ //g' -e 's/\*$//' | sed -e 's/$/\n/' \
                    >>${working}/${group}_per_gene/${genename}.FAA
                
        done

        printf "\n\tDone with ${genename}\n";
        count=$((count+1));    
done

echo "" >${working}/${group}_per_gene/DONE.txt;
printf "\n\tAll Done with extraction\n\n"

