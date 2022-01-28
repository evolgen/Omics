#!/usr/bin/bash

set -e

fasta=$1

module load raxml gcc samtools

dir_name=$(dirname "$fasta") ;
run_name=$(dirname "$fasta" | sed -e 's/\/$//' -e 's/.*\///') ;

cd ${dir_name} ;
sed -i -e '/^>/! s/!/N/g' $fasta ;
samtools faidx $fasta ;
size=$(cut -f2 ${fasta}.fai | head -n 1) ;

echo "" | awk -v size="$size" '{ print "DNA, p1=1-"size"\\3\nDNA, p2=2-"size"\\3\nDNA, p3=3-"size"\\3" }' >cds.partition.txt ;

if [ ! -e "RAxML_bestTree.${run_name}" ]; then
rm -fr RAxML_*${run_name}* ;
printf "\t${run_name} : Started\n" ;
raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRGAMMA -# 100 -p 54321 -q cds.partition.txt -s ${fasta} -n ${run_name} 1>log_raxml 2>&1 ;
printf "\t${run_name} : Done \n" ;
fi

