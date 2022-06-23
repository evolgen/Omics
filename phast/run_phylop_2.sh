#!/usr/bin/bash

file1=$@
# ls file1 in mod_mafs/Lvir_*.maf | parallel -j 28 run_phylop_2.sh

#for file1 in mod_mafs/Lvir_*.maf; do 

file2=$(echo $file1 | sed -e 's/mod_mafs/phylop_scores/' -e 's/\..*/.scores.wig/')
phyloP --wig-scores --method LRT init.mod --msa-format MAF $file1 >$file2


#; done 


