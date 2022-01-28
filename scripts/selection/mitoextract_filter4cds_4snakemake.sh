#!/usr/bin/bash

set -e -o pipefail

module load samtools trimal mafft raxml gcc muscle java macse;
#conda activate hyphy;

infile=$1
genefile=$2

samtools faidx ${infile} ;

inconsistent_count=$(cut -f2 ${infile}.fai | awk '$0%3 != 0' | wc -l);

if [ "$inconsistent_count" -ne "0" ]; then
    printf "\t CDS seq triplet problem : ${infile} = $inconsistent_count sequences\n";
    awk -F'\t' '$2%3 != 0' ${infile}.fai >${infile}.problematic ;
    printf "\t Please check : ${infile}.problematic \n\n" ;
    exit 1;
fi

printf "" >${genefile} ;

grep '^>' $infile | grep -e 'Sebastes_' |
    sed -e 's/^>//' -e 's/ .*//' |
    sort -u |
    while read speciesname; do
        printf ">${speciesname}\n" >> ${genefile} ;
        samtools faidx $infile ${speciesname} |
            grep -v '^>' | 
            sed -e '/^$/d' -e '$s/TAA$/NNN/' -e '$s/TAG$/NNN/' -e '$s/TGA$/NNN/' -e '$s/taa$/NNN/' -e '$s/tag$/NNN/' -e '$s/tga$/NNN/' |
            tr "\n" " " | sed -e 's/ //g' -e 's/$/\n/' |
            fold -w 3 |
            sed -e 's/^TAA$/NNN/' -e '$s/^TAG$/NNN/' -e '$s/^TGA$/NNN/' -e '$s/^AGA$/NNN/' -e '$s/^AGG$/NNN/'
                -e '$s/^taa$/NNN/' -e '$s/^tag$/NNN/' -e '$s/^tga$/NNN/' -e 's/^aga$/NNN/' -e 's/^agg$/NNN/' | 
            awk 'length($0) == 3' |
            tr "\n" " " | sed -e 's/ //g' -e 's/$/\n/' >>${genefile};
    done          

