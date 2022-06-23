#!/usr/bin/bash

set -e 

mkdir -p CHUNKS            # put fragments here
for file in chr*.maf ; do
	root=`basename $file .maf`
	msa_split $file --in-format MAF --refseq $root.fa --windows 1000000,0 --out-root CHUNKS/$root --out-format SS --min-informative 1000 --between-blocks 5000 
done


