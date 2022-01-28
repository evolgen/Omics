#!/usr/bin/bash

set -e

fasta=$1

module load raxml gcc samtools

dir_name=$(dirname "$fasta") ;
run_name=$(dirname "$fasta" | sed -e 's/\/$//' -e 's/.*\///') ;

cd ${dir_name} ;
sed -i -e '/^>/! s/!/X/g' -e '/^>/! s/?/X/g' $fasta ;
samtools faidx $fasta ;
size=$(cut -f2 ${fasta}.fai | head -n 1) ;

echo "" | awk -v size="$size" '{ print "WAG, p1=1-"size" }' >>pep.partition.txt ;

if [ ! -e "RAxML_bestTree.${run_name}" ]; then
rm -fr RAxML_*${run_name}* ;
printf "\t${run_name} : Started\n" ;
raxmlHPC-PTHREADS-SSE3 -T $(eval nproc) -m GTRGAMMA -# 100 -p 54321 -q pep.partition.txt -s ${fasta} -n ${run_name} 1>log_raxml 2>&1 ;
printf "\t${run_name} : Done \n" ;
fi

