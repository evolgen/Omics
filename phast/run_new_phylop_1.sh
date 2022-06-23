#!/usr/bin/bash

file1=$@
# ls file1 in new_Chunks/Lvir_*.ss | parallel -j 28 run_new_phylop_1.sh


file2=$(echo $file1 | sed -e 's/new_Chunks/new_phylop_elements/' -e 's/\..*/.element-scores.txt/')
file3=$(echo $file1 | sed -e 's/new_Chunks/new_gffs/' -e 's/\..*/.gff/')
file4=$(echo $file1 | sed -e 's/new_Chunks/new_phylop_scores/' -e 's/\..*/.scores.wig/')


phyloP --features $file3 --method SCORE --mode CONACC Lacerta_init.mod --msa-format SS $file1 >$file2
phyloP --wig-scores --method LRT Lacerta_init.mod --msa-format MAF $file1 >$file2


#; done 


