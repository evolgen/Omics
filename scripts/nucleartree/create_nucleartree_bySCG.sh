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

if [ ! -e "~/RGP/scripts/nucleartree/list_scp_${group}-genes.filt.txt" ] ; then

    printf " Creating the list of genomes\n";

    for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/; do 
        name=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/\/.*//'); 
        printf $name"\t"; 
        ls $file1/*.fna | wc -l; 
    done >~/RGP/scripts/nucleartree/list_${group}_genomes.txt

    awk -F'\t' -v limit="${limit}" '$2>=limit {print $1}' ~/RGP/scripts/nucleartree/list_${group}_genomes.txt |
        sed -e '/^$/d' \
        >~/RGP/scripts/nucleartree/list_nucleargenomes_${group}.txt
    awk -F'\t' '$2<150 {print $1}' ~/RGP/scripts/nucleartree/list_${group}_genomes.txt |
        sed -e '/^$/d' \
        >~/RGP/scripts/nucleartree/filtered_nucleargenomes_${group}.txt

    number=$(wc -l ~/RGP/scripts/nucleartree/list_nucleargenomes_${group}.txt | awk '{print $1}')

    printf " Creating the list of genes for ${number} genomes\n";
    \ls /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/*.fna |
        fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes_${group}.txt |
        sed -e 's/\/.*\///' |
        sort | uniq -c |
        awk -v number="$number" '$1==number' |
        sed -e '/^$/d' \
            >~/RGP/scripts/nucleartree/list_scp_${group}-genes.txt

    awk '{print $2}' ~/RGP/scripts/nucleartree/list_scp_${group}-genes.txt \
        >~/RGP/scripts/nucleartree/list_scp_${group}-genes.filt.txt

fi
    

printf " Starting the concatenation of genes for each species\n";

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree"
mkdir -p ${working}/${group}_per_species

count=1
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/ |
    fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes_${group}.txt |
    while read folder; do
        speciesname=$(echo ${folder} | sed -e 's/.*\/output\///' -e 's/\/.*//');
        printf " Running for ${speciesname} at\t${folder}\n\t";
        printf ">${speciesname}\n" \
            >${working}/${group}_per_species/${speciesname}.FAS

        cat ~/RGP/scripts/nucleartree/list_scp_${group}-genes.filt.txt | 
            while read gene; do
                fastaname=$(echo ${gene} | sed -e 's/.*\///' -e 's/\.fna//');    
                printf "${fastaname}  ";

                if [ ! -e "${folder}/${gene}" ]; then
                    printf "ERROR : Missing in ${speciesname}\n ${folder}/${gene}\n\n";
                    exit 1
                fi
                    
                cat ${folder}/${gene} |
                    grep -v '^>' |
                    tr "\n" " " | 
                    sed -e 's/ //g' |
                    sed -e 's/taa$//' -e 's/tag$//' -e 's/tga$//' >>${working}/${group}_per_species/${speciesname}.FAS 

                done

                printf "\n" >>${working}/${group}_per_species/${speciesname}.FAS
                printf "\nDone with ${speciesname}\n";

            done


\ls -d ${working}/${group}_per_species/*.FAS |
    xargs cat |
    sed -e 's/>/\n>/' |
    sed -e '/^$/d' \
        >${working}/${group}_CDS.FASTA
sed -i -e '/^$/d' ${working}/${group}_CDS.FASTA;


module load trimal mafft raxml gcc muscle java macse prank

printf "\n\tUsing PRANK for MSA\n";
prank -verbose -d=${working}/${group}_CDS.FASTA -o=${working}/${group}_CDS.prank.FASTA -DNA -F 

printf "\tUsing TRIMAL for trimming the MSA\n";
trimal -in ${working}/${group}_CDS.prank.FASTA.best.fas -fasta -out ${working}/${group}_CDS.pranktrim.FASTA -gt 1 


cd ${working}
printf "\n\tUsing RAxML for tree building\n";
raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 01859 -x 01859 -# 100 -s ${working}/${group}_CDS.pranktrim.FASTA -n Rockfish_${group}_nuclear -o Adelosebastes_latens



