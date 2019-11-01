#!/usr/bin/bash

set -e

# ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/alignments/wga/masurca/*/align*.paf | parallel -j 32 sh /global/home/users/rohitkolora/RGP/snakemake_workflows/illumina/align_genomes/parallel_masurca_divscore_4m_pafs.sh

file1=$@

    species1=$(echo $file1 | sed -e 's/.*\///' -e 's/align\.//' -e 's/\..*//');
	species2=$(echo $file1 | sed -e 's/.*\///' -e 's/align\.//' -e 's/\.paf//' -e 's/.*\.//');
    name=$(basename $file1)
	printf $species1"\t"$species2"\t" >${file1/%.paf/.score};
##	awk -F'\t' '$NF>=12 && $12==60 && $11>0 {score+=($10/$11); count+=1} END {print score/count}' $file1 >>${file1/%.paf/.score};

    python3 ~/RGP/snakemake_workflows/illumina/align_genomes/divergence_frompaf.py ${file1} >${file1/%.paf/.identity};

    awk '{ mismatch+=$(NF-4); total+=$(NF-2) } END { print mismatch"\t"total }' ${file1/%.paf/.identity} | 
        awk -F'\t' '{print $2/$1}' \
        >>${file1/%.paf/.score};

    printf "Done with ${name}\n";
