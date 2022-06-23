#!/usr/bin/bash

set -e 

: '
mkdir -p CHUNKS            # put fragments here
for file in chr*.maf ; do
	root=`basename $file .maf`
	msa_split $file --in-format MAF --refseq $root.fa --windows 1000000,0 --out-root CHUNKS/$root --out-format SS --min-informative 1000 --between-blocks 5000 
done
'

file=$@
mkdir -p TREES     # put estimated tree models here
rm -f TREES/*      # in case old versions left over
#for file in chr*.*.ss ; do 
	root=`basename $file .ss` 
	phastCons --target-coverage 0.125 --expected-length 20 --gc 0.4 --estimate-trees TREES/$root $file init.mod --no-post-probs
done

