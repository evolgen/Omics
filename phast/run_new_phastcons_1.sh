#!/usr/bin/bash

mkdir -p new_Chunks new_Trees            # put fragments here

file=$@

#for file in mod_mafs/*.maf
#do
	root1=$(echo $file | sed -e 's/new_mod_mafs/new_Chunks/' -e 's/\.maf//')
	seq=$(echo $file | sed -e 's/new_mod_mafs/seqs/' -e 's/\.maf/.fa/')
	msa_split $file --in-format MAF --refseq $seq --windows 5000000,0 --out-root $root1 --out-format SS --min-informative 30 --between-blocks 1000 
	echo "$root	DONE"
	
#done


