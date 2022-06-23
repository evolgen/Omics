#!/usr/bin/bash

file1=$@
# ls file1 in new_mod_mafs/Lvir_*.maf | parallel -j 28 run_new_phylop_2.sh


file2=$(echo $file1 | sed -e 's/new_mod_mafs/new_phylop_scores/' -e 's/\..*/.scores.wig/')
phyloP --wig-scores --method LRT Lacerta_init.mod --msa-format MAF $file1 >$file2


#; done 


