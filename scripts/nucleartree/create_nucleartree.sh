#!/usr/bin/bash

set -e -o pipefail

module load seqtk samtools bedtools

if [ ! -e "~/RGP/scripts/nucleartree/list_scp_cvg-genes.filt.txt" ] ; then

    printf " Creating the list of genomes\n";
    awk -F'\t' '$2>=150 {print $1}' /global/scratch2/rohitkolora/Rockfish/Genomes/estimates/cvg_genomes.txt |
        sed -e '/^$/d' \
        >~/RGP/scripts/nucleartree/list_nucleargenomes.txt
    awk -F'\t' '$2<150 {print $1}' /global/scratch2/rohitkolora/Rockfish/Genomes/estimates/cvg_genomes.txt |
        sed -e '/^$/d' \
        >~/RGP/scripts/nucleartree/filtered_nucleargenomes.txt

    number=$(wc -l ~/RGP/scripts/nucleartree/list_nucleargenomes.txt | awk '{print $1}')

    printf " Creating the list of genes for ${number} genomes\n";
    \ls /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_cvg/single_copy_busco_sequences/*.fna |
        fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes.txt |
        sed -e 's/\/.*\///' |
        sort | uniq -c |
        awk -v number="$number" '$1==number' |
        sed -e '/^$/d' \
            >~/RGP/scripts/nucleartree/list_scp_cvg-genes.txt

    awk '{print $2}' ~/RGP/scripts/nucleartree/list_scp_cvg-genes.txt \
        >~/RGP/scripts/nucleartree/list_scp_cvg-genes.filt.txt

fi
    

printf " Starting the concatenation of genes for each species\n";

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree"
mkdir -p ${working}/per_species

count=1
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_cvg/single_copy_busco_sequences/ |
    fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes.txt |
    while read folder; do
        speciesname=$(echo ${folder} | sed -e 's/.*\/output\///' -e 's/\/.*//');
        printf " Running for ${speciesname} at\t${folder}\n\t";
        printf ">${speciesname}\n" \
            >${working}/per_species/${speciesname}.FAS

        cat ~/RGP/scripts/nucleartree/list_scp_cvg-genes.filt.txt | 
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
                    sed -e 's/taa$//' -e 's/tag$//' -e 's/tga$//' >>${working}/per_species/${speciesname}.FAS 

                done

                printf "\n" >>${working}/per_species/${speciesname}.FAS
                printf "\nDone with ${speciesname}\n";

            done


\ls -d ${working}/per_species/*.FAS |
    xargs cat |
    sed -e 's/>/\n>/' |
    sed -e '/^$/d' \
        >${working}/CDS.FASTA
sed -i -e '/^$/d' ${working}/CDS.FASTA;


module load trimal mafft raxml gcc muscle java macse prank

printf "\n\tUsing PRANK for MSA\n";
prank -verbose -d=${working}/CDS.FASTA -o=${working}/CDS.prank.FASTA -DNA -F 

printf "\tUsing TRIMAL for trimming the MSA\n";
trimal -in ${working}/CDS.prank.FASTA.best.fas -fasta -out ${working}/CDS.pranktrim.FASTA -gt 1 


cd ${working}
printf "\n\tUsing RAxML for tree building\n";
raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 01859 -x 01859 -# 100 -s ${working}/CDS.pranktrim.FASTA -n Rockfish_nuclear -o Adelosebastes_latens



