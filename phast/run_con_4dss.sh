#!/bin/bash

set -e

file1=$@

# for file1 in ./mafs/Lvir_*.maf; do 
infile=$(echo $file1 | sed -e 's/mod_mafs/seqs/' -e 's/\.maf/.fa/')
file2=$(echo $file1 | sed -e 's/mod_mafs/gffs/' -e 's/\.maf/.gff/')
file3=$(echo $file2 | sed -e 's/gffs/codons_ss/' -e 's/\.gff/.codons.ss/')
file4=$(echo $file3 | sed -e 's/codons_ss/sites_ss/' -e 's/\.ss/.sites.ss/')


if [ -f $file2 ]
then

msa_view $file1 --in-format MAF --features $file2 --catmap "NCATS = 3; CDS 1-3" --out-format SS --unordered-ss --reverse-groups gene_id >$file3
msa_view $file3 --in-format SS --out-format SS --tuple-size 1 > $file4

fi


