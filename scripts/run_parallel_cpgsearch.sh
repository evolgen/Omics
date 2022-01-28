#!/usr/bin/bash

set -e

fasta=$@

module load emboss
echo $fasta;
mkdir -p $(dirname "$fasta")/CpG && cd $(dirname "$fasta")/CpG ;
filename=$(basename "$fasta" | sed -e 's/.*\///')
#cpgreport -sequence $fasta -score 17 -outfile ${filename}.cpgreport.txt -outfeat ${filename}.GC.gff
newcpgreport -sequence $fasta -outfile ${filename}.cpgreport.txt -window 200 -shift 1 -minlen 200 -minoe 0.6 -minpc 0.5


