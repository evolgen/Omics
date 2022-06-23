#!/usr/bin/bash

file1=$@
# ls file1 in new_mod_mafs/Lvir_*.maf | parallel -j 28 run_new_phylop_maf.sh

#rm new_phylop_elements/* new_phylop_scores/* new_phylop_elements_bb/*

file2=$(echo $file1 | sed -e 's/new_mod_mafs/new_phylop_elements_other/' -e 's/\..*/.element-scores.txt/')
file3=$(echo $file1 | sed -e 's/new_mod_mafs/other_gffs/' -e 's/\..*/.gff/')


phyloP --features $file3 --method SCORE --mode CONACC Lacerta_init.mod --msa-format MAF $file1 >$file2



#; done 


