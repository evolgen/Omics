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
prank -verbose -d=${working}/${group}_CDS.FASTA -o=${working}/${group}_CDS.prank.FASTA -DNA -F

printf "\tUsing TRIMAL for trimming the MSA\n";
trimal -in ${working}/${group}_CDS.prank.FASTA.best.fas -fasta -out ${working}/${group}_CDS.pranktrim.FASTA -gt 1

cd ${working}
printf "\n\tUsing RAxML for tree building\n";
raxmlHPC-PTHREADS-SSE3 -T 24 -f a -m GTRGAMMA -p 185 -x 185 -# 100 -s ${working}/${group}_CDS.pranktrim.FASTA -n Rockfish_${group}_pranknuclear -o Sebastolobus_alascanus,Sebastolobus_altivelis


