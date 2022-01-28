#!/usr/bin/bash

set -e

if [ "$#" -ne "3" ]; then
    echo "Need inputs : Reference Query Name"
    exit 1;
fi

module load minimap2 

reference=$1
query=$2
name=$3

minimap2 -x asm20 -t $(eval nproc) --cs=long -Lc $reference $query >${name}.paf
python3 ~/RGP/snakemake_workflows/illumina/align_genomes/divergence_frompaf.py ${name}.paf >${name}.identity
awk '{ mismatch+=$(NF-4); total+=$(NF-2) } END { print mismatch"\t"total }' ${name}.identity | 
    awk -F'\t' '{print $2/$1}' >${name}.score

printf "${name}\t${reference}\t${query}\t"; cat ${name}.score;
echo

