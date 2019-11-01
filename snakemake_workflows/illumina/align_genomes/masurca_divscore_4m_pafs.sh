#!/usr/bin/bash

set -e

for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*/align*.paf; do

	species1=$(echo $file1 | sed -e 's/.*\///' -e 's/align\.//' -e 's/\..*//');
	species2=$(echo $file1 | sed -e 's/.*\///' -e 's/align\.//' -e 's/\.paf//' -e 's/.*\.//');
	printf $species1"\t"$species2"\t" >${file1/%.paf/.score};
##	awk -F'\t' '$NF>=12 && $12==60 && $11>0 {score+=($10/$11); count+=1} END {print score/count}' $file1 >>${file1/%.paf/.score};

    python3 ~/RGP/snakemake_workflows/illumina/align_genomes/divergence_frompaf.py ${file1} >${file1/%.paf/.identity};
    awk -F'\t' '{sum+=$NF; count+=1} END {print sum/count}' $file1 >>${file1/%.paf/.score};

done	

for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*_*; do

	species=$(echo $file1 | sed -e 's/.*\///')
	printf $species"\t"$species"\t"1"\n" >${file1}/align.${species}.${species}.score;

done	

cat /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*/*.score >/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/stats/divergence_masurca_wga.txt;
cat /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*/*.score | awk -F'\t' 'BEGIN{OFS="\t"} {$NF=(1-$NF)*100; print $0}' >/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/stats/dissimilarity_masurca_wga.txt;


