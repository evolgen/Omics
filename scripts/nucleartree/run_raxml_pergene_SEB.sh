#!/usr/bin/sh

set -e -o pipefail

file1=$@
working=$(dirname $file1)

name=$(echo $file1 | sed -e 's/.*\///' -e 's/\.FAS$//')
file_prank=$(echo $file1 | sed -e 's/\.FASTA//' -e 's/\.fasta//' -e 's/\.fas$//' -e 's/\.FAS$//' -e 's/\.fa$//' -e 's/\.FA$//' -e 's/.prank.fas//')
file_trim=$(echo $file_prank | sed -e 's/$/.pranktrim.FASTA/')

module load trimal mafft raxml gcc muscle java macse prank

if [ ! -e "${file_prank}" ]; then
    printf "\n\tUsing PRANK for name\n";
    prank -verbose -d=${file1} -o=${file_prank} -DNA -F
fi

if [ ! -e "${file_trim}" ]; then
    printf "\tUsing TRIMAL for trimming the MSA\n";
    trimal -in ${file_prank}.best.fas -fasta -out ${file_trim} -gt 1
fi

cd ${working}
number=$(shuf -i 1-10000 -n 1)
printf "\n\tUsing RAxML for tree building \t $name\n";
raxmlHPC-SSE3 -T 5 -f a -m GTRGAMMA -p $number -x $number -# 5 -s ${file_trim} -n RF_acti_SCG_${name} -o Sebastolobus_alascanus


