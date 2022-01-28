#!/usr/bin/bash

set -e -o pipefail

find  /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/large_sequence_files/ -type f -name "seq.FNA.fai" | 
    xargs fgrep -c -f ~/RGP/scripts/selection/sebastes_species.alascanus.list | 
    sed -e 's/:/\t/' | 
    awk '$2==77 {gsub(".fai","",$1); print $1}' | 
    parallel -j 32 bash ~/RGP/scripts/trees/parallel_lifted_77species_copy.sh

