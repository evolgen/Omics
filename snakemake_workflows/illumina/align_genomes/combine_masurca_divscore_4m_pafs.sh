#!/usr/bin/sh

set -e 

for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*_*; do
        species=$(echo $file1 | sed -e 's/.*\///')
        printf $species"\t"$species"\t0\n" >${file1}/align.${species}.${species}.score;
done

cat /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*/*.score >/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/stats/divergence_masurca_wga.txt;
cat /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*/*.score | awk -F'\t' '{print $1"\t"$2"\t"$3*100}' >/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/stats/dissimilarity_masurca_wga.txt;


