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

number=$(wc -l ~/RGP/scripts/nucleartree/Latest/list_nucleargenomes_${group}.uniq.latest.txt | awk '{print $1}');
printf " Creating the list of genes for ${number} unique genomes\n";

scp_genes=$(cat ~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt | wc -l);

printf "\tSCO genes in total\t${scp_genes}\n";
minimal=$(($limit/10));

if [ "$scp_genes" -lt "${minimal}" ]; then
    printf "\n\tToo Few SCO genes - ${scp_genes}\n\n";
    ##rm -f ~/RGP/scripts/nucleartree/Latest/*${group}*.latest.txt;
    exit 0;
fi    

printf " Starting the concatenation of genes for each species\n";

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/BUSCO"
mkdir -p ${working}/${group}_per_species;


if [ ! -e "${working}/${group}_PEP.FASTA" ]; then
count=1
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/ |
    fgrep -w -f ~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt |
    while read folder; do
        speciesname=$(echo ${folder} | sed -e 's/.*\/output\///' -e 's/\/.*//');
        printf " Running for ${speciesname} at\t${folder}\n\t";
        printf ">${speciesname}\n" \
            >${working}/${group}_per_species/${speciesname}.PEP

        cat ~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt |
            sed -e 's/\.fna$/.faa/' | 
            while read gene; do
                fastaname=$(echo ${gene} | sed -e 's/.*\///' -e 's/\.faa//');    
                printf "${fastaname}  ";

                if [ ! -e "${folder}/${gene}" ]; then
                    printf "ERROR : Missing in ${speciesname}\n ${folder}/${gene}\n\n";
                    exit 1
                fi
                    
                cat ${folder}/${gene} |
                    grep -v '^>' |
                    tr "\n" " " | 
                    sed -e 's/ //g' |
                    sed -e 's/\*$//' -e 's/\*/X/g' \
                        >>${working}/${group}_per_species/${speciesname}.PEP 
                done

                printf "\n" >>${working}/${group}_per_species/${speciesname}.PEP;
                printf "\nDone with ${speciesname}\n";

            done
fi


\ls -d ${working}/${group}_per_species/*.PEP |
    xargs cat |
    sed -e 's/>/\n>/' |
    sed -e '/^$/d' \
        >${working}/${group}_PEP.FASTA
sed -i -e '/^$/d' ${working}/${group}_PEP.FASTA;


module load trimal mafft raxml gcc muscle java 

printf "\n\tUsing MAFFT for MSA\n";
mafft --auto ${working}/${group}_PEP.FASTA >${working}/${group}_PEP.mafft.FASTA 

printf "\tUsing TRIMAL for trimming the MSA\n";
trimal -in ${working}/${group}_PEP.mafft.FASTA -fasta -out ${working}/${group}_PEP.maffttrim.FASTA -gt 1 

cd ${working}
printf "\n\tUsing RAxML for tree building\n";
raxmlHPC-PTHREADS-SSE3 -T 32 -f a -m PROTGAMMAGTR -p 01859 -x 01859 -# 10000 -s ${working}/${group}_PEP.maffttrim.FASTA -n Rockfish_${group}_pep -o Sebastolobus_alascanus,Sebastolobus_altivelis



