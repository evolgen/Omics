#!/usr/bin/bash

set -e -o pipefail

if [ "$#" -ne 3 ]; then
    echo "    Need 3 input arguments";
    exit 1;
fi

log_file=$1
qry_gen=$2
qry_gen_rc=$3

printf "" >${qry_gen_rc} ;
cat ${log_file} | while read num; do
    chromo=$(awk -v num="$num" 'NR==num' ${qry_gen}.fai | cut -f1);
    printf "\tReverse complementing  :  ${chromo}\t${num}\n";
    printf ">${num}\n" >>${qry_gen_rc} ;                 # RC header derived
    samtools faidx ${qry_gen} ${num} |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${qry_gen_rc} ;            # RC converted
    done    

