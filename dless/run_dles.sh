#!/usr/bin/bash

#ls ../mod_mafs/Lvir_*.maf

file1=$@

#for file1 in ../mod_mafs/Lvir_*.maf; do 

name=$(echo $file1 | sed -e 's/\.maf/.gff/' -e 's/.*Lvir_/Lvir_/')
file2=$(echo $file1 | sed -e 's/\.\.\/mod_mafs/.\/histories/' -e 's/\.maf/.ih/')
file3=$(echo $file1 | sed 's/\.\.\/mod_mafs\//.\/models\//' | sed 's/.maf/.model/')
echo "dless -H $file2 --msa-format MAF $file1 $file3 >./dless/$name"
dless -H $file2 --msa-format MAF $file1 $file3 >./dless/$name

#;done


