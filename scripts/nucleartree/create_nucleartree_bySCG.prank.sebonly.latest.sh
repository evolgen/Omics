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

working="/global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/BUSCO"

module load trimal prank raxml gcc muscle java 
printf "\n\tUsing PRANK for MSA\n";
prank -verbose -d=${working}/sebonly_${group}_CDS.FASTA -o=${working}/sebonly_${group}_CDS.prank.FASTA -DNA -F

printf "\tUsing TRIMAL for trimming the MSA\n";
trimal -in ${working}/sebonly_${group}_CDS.prank.FASTA.best.fas -fasta -out ${working}/sebonly_${group}_CDS.pranktrim.FASTA -gt 1

cd ${working}
printf "\n\tUsing RAxML for tree building\n";
raxmlHPC-PTHREADS-SSE3 -T 32 -f a -m GTRGAMMA -p 189 -x 189 -# 100 -s ${working}/sebonly_${group}_CDS.pranktrim.FASTA -n Sebastesonly_${group}_pranknuclear -o Sebastiscus_albofasciatus


