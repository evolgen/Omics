#!/usr/bin/bash

file1=$@
# ls file1 in Chunks/Lvir_*.ss | parallel -j 28 run_phylop_1.sh

#for file1 in Chunks/Lvir_*.ss; do 

file2=$(echo $file1 | sed -e 's/Chunks/phylop_elements/' -e 's/\..*/.element-scores.txt/')
file3=$(echo $file1 | sed -e 's/Chunks/gffs/' -e 's/\..*/.gff/')
phyloP --features $file3 --method SCORE --mode CONACC init.mod --msa-format SS $file1 >$file2

#; done 


