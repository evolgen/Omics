#!/usr/bin/bash

#ls ../mod_mafs/Lvir_*.maf

file=$@

#for file1 in ../mod_mafs/Lvir_*.maf; do 

file2=$(echo $file | sed -e 's/.ss/.ih/' -e 's/.*Lvir_/Lvir_/' -e 's/\.1-[0-9]*\././')
file3=$(echo $file2 | sed 's/\.ih$/.model/')
file4=$(echo $file2 | sed 's/\.ih$/.tree/')
cat ./models/$file3 | grep TREE | sed 's/TREE: //' >./trees/$file4
indelHistory --msa-format SS $file ./trees/$file4 >./histories/${file2}



