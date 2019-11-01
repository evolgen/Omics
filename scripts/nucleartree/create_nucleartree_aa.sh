#!/usr/bin/bash

set -e -o pipefail

module load seqtk samtools bedtools

sed -e 's/\.fna/.faa/' ~/RGP/scripts/nucleartree/list_scp_cvg-genes.filt.txt >~/RGP/scripts/nucleartree/list_scp_cvg-genes.filt.aa.txt
printf " Starting the concatenation of genes for each species\n";

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/AA"
mkdir -p ${working}/per_species_AA


count=1
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/BUSCO/run_busco_cvg/single_copy_busco_sequences/ |
    fgrep -f ~/RGP/scripts/nucleartree/list_nucleargenomes.txt |
    while read folder; do
        speciesname=$(echo ${folder} | sed -e 's/.*\/output\///' -e 's/\/.*//');
        printf " Running for ${speciesname} at\t${folder}\n\t";
        printf ">${speciesname}\n" \
            >${working}/per_species_AA/${speciesname}.FAS

        cat ~/RGP/scripts/nucleartree/list_scp_cvg-genes.filt.aa.txt | 
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
                    sed -e 's/taa$//' -e 's/tag$//' -e 's/tga$//' >>${working}/per_species_AA/${speciesname}.FAS 

                done

                printf "\n" >>${working}/per_species_AA/${speciesname}.FAS
                printf "\nDone with ${speciesname}\n";

            done


\ls -d ${working}/per_species_AA/*.FAS |
    xargs cat |
    sed -e 's/>/\n>/' |
    sed -e '/^$/d' \
        >${working}/CDS.AA.FASTA
sed -i -e '/^$/d' ${working}/CDS.AA.FASTA;


module load trimal mafft raxml gcc muscle java macse prank

printf "\n\tUsing PRANK for MSA\n";
prank -verbose -d=${working}/CDS.AA.FASTA -o=${working}/CDS.AA.prank.FASTA -protein -F 

printf "\tUsing TRIMAL for trimming the MSA\n";
trimal -in ${working}/CDS.AA.prank.FASTA.best.fas -fasta -out ${working}/CDS.AA.pranktrim.FASTA -gt 1 


cd ${working}
printf "\n\tUsing RAxML for tree building\n";
raxmlHPC-PTHREADS-SSE3 -T 32 -f a -m PROTGAMMAAUTO -p 122 -x 122 -# 10000 -s ${working}/CDS.AA.pranktrim.FASTA -n Rockfish_nuclear_AA10k -o Adelosebastes_latens



