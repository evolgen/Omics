#!/usr/bin/bash

set -e -o pipefail

# find  /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/large_sequence_files/ -type f -name "seq.FNA" | parallel -j 32 bash ~/RGP/scripts/trees/obtain_77species_files.large.select.sh

module load samtools 

file1=$@

outdir=$(dirname "$file1" | sed -e 's/\/large_sequence_files/\/Sebastes_77\/select_sequence_files/')
mkdir -p $outdir

cat ${file1}.fai | 
    grep -f <(cat ~/RGP/scripts/selection/sebastes_species.alascanus.list | sed -e 's/^/^/') | 
    cut -f1 |   
    xargs samtools faidx $file1 >${outdir}/seq.FNA ;

echo $outdir


