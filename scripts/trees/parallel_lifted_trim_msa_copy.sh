#!/usr/bin/bash

set -e -o pipefail

pep=$@

module load samtools

if [ "$#" -ne 1 ]; then
    echo "Need atleast one argument"
    exit 1;
fi

#while read pep; do 

cds=$(echo $pep | sed -e 's/pep.seq.msa.trim/cds.seq.msa.trim/'); 
name=$(dirname "$pep" | sed -e 's/.*\///'); 

echo $name; 
##sed -e '/^>/! s/!/X/g' -e '/^>/! s/?/X/g' -e '/^>/ s/.FUN*//' $pep >all_seqs/$name.faa; 
##sed -e '/^>/! s/!/X/g' -e '/^>/! s/?//g' -e '/^>/ s/.FUN*//' $cds >all_seqs/$name.fna; 

##samtools faidx all_seqs/$name.fna; 
##samtools faidx all_seqs/$name.faa; 

cat ~/RGP/scripts/selection/sebastes_species.alascanus.list |
    xargs samtools faidx all_seqs/$name.fna >all_seqs_filt/$name.fna ;
samtools faidx all_seqs_filt/$name.fna ;    
cat ~/RGP/scripts/selection/sebastes_species.alascanus.list |
    xargs samtools faidx all_seqs/$name.faa >all_seqs_filt/$name.faa ;
samtools faidx all_seqs_filt/$name.faa ;    

#done | head -n 2


