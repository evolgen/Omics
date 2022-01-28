#!/usr/bin/bash

set -e -o pipefail

module load samtools trimal mafft raxml gcc muscle java macse emboss;
#conda activate hyphy;

infile=$1
genefile=$2

samtools faidx ${infile} ;

inconsistent_count=$(cut -f2 ${infile}.fai | awk '$0%3 != 0' | wc -l);

if [ "$inconsistent_count" -ne "0" ]; then
    printf "\t CDS seq triplet problem : ${infile} = $inconsistent_count sequences\n";
    awk -F'\t' '$2%3 != 0' ${infile}.fai >${infile}.problematic ;
    printf "\t Please check : ${infile}.problematic \n\n" ;
fi

printf "" >${genefile} ;

grep '^>' $infile | #grep -e 'Sebastes_' -e 'Sebastolobus_alascanus' |
    sed -e 's/^>//' -e 's/ .*//' |
    sort -u |
    while read speciesname; do
        printf ">${speciesname}\n" >> ${genefile} ;
        samtools faidx $infile ${speciesname} |
            grep -v '^>' >>${genefile}; 
    done          

/global/scratch2/rohitkolora/miniconda3/envs/hyphy/bin/hmmer2go getorf -i ${genefile} -o ${genefile}.orfout -c -t 2 --nomet ;
sed -e '/^>/ s/_[0-9]* \[[0-9]* - [0-9]*\] (.*//' ${genefile}.orfout >${genefile}.orf



