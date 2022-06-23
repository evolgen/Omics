#!/usr/bin/bash

mkdir -p Chunks Trees            # put fragments here

file=$@

#for file in mod_mafs/*.maf
#do
	root1=$(echo $file | sed -e 's/mod_mafs/Chunks/' -e 's/\.maf//')
	seq=$(echo $file | sed -e 's/mod_mafs/seqs/' -e 's/\.maf/.fa/')
	msa_split $file --in-format MAF --refseq $seq --windows 5000000,0 --out-root $root1 --out-format SS --min-informative 1000 --between-blocks 5000 
	echo "$root	DONE"
	
#done



