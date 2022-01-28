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

###if [ ! -e "~/RGP/scripts/nucleartree/Latest/list_${group}_genomes.latest.txt" ]; then
###    printf " Creating the list of genomes\n";
###    \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/ |
###    grep -v 'SEB-74' |
###    while read file1; do
###        name=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/\/.*//');
###        sample=$(echo $file1 | sed -e 's/.*\/output\///' -e 's/\/assembly\/.*//' -e 's/.*\///');
###        printf $name"\t"$sample"\t";
###        \ls -d $file1/*.fna | wc -l;
###    done >~/RGP/scripts/nucleartree/Latest/list_${group}_genomes.latest.txt
###fi

if [ ! -e "~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt" ] ; then

    awk -F'\t' -v limit="${limit}" '$3>=limit {print $1}' ~/RGP/scripts/nucleartree/Latest/list_${group}_genomes.latest.txt |
        sed -e '/^$/d' \
        >~/RGP/scripts/nucleartree/Latest/list_nucleargenomes_${group}.latest.txt;
    awk -F'\t' -v limit="${limit}" '$3<limit {print $1}' ~/RGP/scripts/nucleartree/Latest/list_${group}_genomes.latest.txt |
        sed -e '/^$/d' \
        >~/RGP/scripts/nucleartree/Latest/filtered_nucleargenomes_${group}.latest.txt;

    sort -k1,1 -k2,2 -k3,3nr ~/RGP/scripts/nucleartree/Latest/list_${group}_genomes.latest.txt |
        awk -F'\t' '!x[$1]++ {print $1"\t"$2"\t"$3}' |
        fgrep -f ~/RGP/scripts/nucleartree/selected_species.txt \
        >~/RGP/scripts/nucleartree/Latest/list_nucleargenomes_${group}.uniq.latest.txt;
    awk -F'\t' '{print $2}' ~/RGP/scripts/nucleartree/Latest/list_nucleargenomes_${group}.uniq.latest.txt \
        >~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt;
fi

number=$(wc -l ~/RGP/scripts/nucleartree/Latest/list_nucleargenomes_${group}.uniq.latest.txt | awk '{print $1}');
printf " Creating the list of genes for ${number} unique genomes\n";

if [ ! -e "~/RGP/scripts/nucleartree/Latest/list_all_samples_${group}.files.latest.lst" ]; then

    find /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_acti/single_copy_busco_sequences/ -type f -name '*.fna' \
        >~/RGP/scripts/nucleartree/Latest/list_all_samples_${group}.files.latest.lst;

fi    

if [ ! -e "~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt" ] ; then

    cat ~/RGP/scripts/nucleartree/Latest/list_all_samples_${group}.files.latest.lst |
        fgrep -f ~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt |
        sed -e 's/\/.*\///' |
        sort | uniq -c |
        awk -v number="$number" '$1==number' |
        sed -e '/^$/d' \
            >~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.latest.txt;

    awk '{print $2}' ~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.latest.txt \
        >~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt;

fi

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

count=1
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_${group}/single_copy_busco_sequences/ |
    fgrep -w -f ~/RGP/scripts/nucleartree/Latest/list_nuclear_samples_${group}.uniq.latest.txt |
    while read folder; do
        speciesname=$(echo ${folder} | sed -e 's/.*\/output\///' -e 's/\/.*//');
        printf " Running for ${speciesname} at\t${folder}\n\t";
        printf ">${speciesname}\n" \
            >${working}/${group}_per_species/${speciesname}.FAS

        cat ~/RGP/scripts/nucleartree/Latest/list_scp_${group}-genes.filt.latest.txt | 
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
                    sed -e 's/taa$//' -e 's/tag$//' -e 's/tga$//' \
                        >>${working}/${group}_per_species/${speciesname}.FAS 

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


module load trimal mafft raxml gcc muscle java 

printf "\n\tUsing MAFFT for MSA\n";
mafft --auto --anysymbol --thread 32 ${working}/${group}_CDS.FASTA >${working}/${group}_CDS.prank.FASTA 

printf "\tUsing TRIMAL for trimming the MSA\n";
trimal -in ${working}/${group}_CDS.prank.FASTA -fasta -out ${working}/${group}_CDS.maffttrim.FASTA -gt 1 


cd ${working}
printf "\n\tUsing RAxML for tree building\n";
raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 01859 -x 01859 -# 10000 -s ${working}/${group}_CDS.maffttrim.FASTA -n Rockfish_${group}_nuclear -o Sebastolobus_alascanus,Sebastolobus_altivelis



