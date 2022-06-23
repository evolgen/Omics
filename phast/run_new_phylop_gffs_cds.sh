#!/usr/bin/bash

set -e

file1=$@
# ls file1 in /scr/k61san/nowicklab/Lacerta/Multiz/new_gffs_cds/Lvir_*.gff | parallel -j 28 sh /scr/k61san/nowicklab/Lacerta/Multiz/run_new_phylop_gffs_cds.sh
# mkdir -p /scr/k61san/nowicklab/Lacerta/Multiz/new_phylop_gffs_cds

file2=$(echo $file1 | sed -e 's/new_gffs_cds/new_phylop_gffs_cds/' -e 's/\.gff$/.score.gff/')
file3=$(echo $file1 | sed -e 's/new_gffs_cds/new_mod_mafs/' -e 's/\.gff$/.maf/')
file4=$(echo $file1 | sed -e 's/new_gffs_cds/new_Chunks/' -e 's/\.gff$/.[0-9]*-[0-9]*.ss/')
file5=$(ls $file4)

#echo $file5

echo "phyloP --features $file1 --method SCORE --mode CONACC Lacerta_init.mod --msa-format SS $file5 >$file2"
/scr/bloodymary/rohit/phast/bin/phyloP --features $file1 --method SCORE --mode CONACC /scr/k61san/nowicklab/Lacerta/Multiz/Lacerta_init.mod --msa-format SS $file5 >$file2


#; done 


