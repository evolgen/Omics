#!/usr/bin/bash

set -e

for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/unordered/PS001_fishes/*_R1_*.fastq.gz; do 
    
    file2=$(echo $file1 | sed -e 's/_R1_/_R2_/')
    pattern=$(echo $file1 | sed -e 's/.*\///' -e 's/_.*//' -e 's/Seb/SEB-/');
    awk -F'\t' -v pattern="$pattern" '$1==pattern' /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/Master_table.tsv | while read line; do
      species=$(echo $line | sed -e 's/ /\t/' | sed -e 's/.*\t//' -e 's/ /_/');
      folder=$(echo $species | sed -e 's/Sebastes_/S-/' -e 's/Sebastolobus_/B-/' -e 's/Helicolenus_/H-/' -e 's/Adelosebastes_/A-/')
      sample=$(echo $line | sed -e 's/ /\t/' | sed -e 's/\t.*//');
      mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/${species}
      mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/${species}/${folder}_$sample
      cp $file1 $file2 /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/${species}/${folder}_$sample
      printf "$species\t\t$sample\n"
    done

done


