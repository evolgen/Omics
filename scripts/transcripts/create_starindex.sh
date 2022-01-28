#!/usr/bin/bash

set -e 

#fasta1=$1

#module load star ;
# conda activate funannotate

#\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/referencegenome.fasta | grep -v -e ruber | while read fasta1; do
\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/Sebastes_ruberrimus/S-ruberrimus_Seb91/assembly/masurca/final.genome.scf.FAS | 
    while read fasta1; do
    echo $fasta1 ;
    cd $(dirname "$fasta1");
    mkdir -p star_index ;
    /global/scratch2/rohitkolora/miniconda3/envs/funannotate/bin/STAR --limitGenomeGenerateRAM 300000000000 --runThreadN $(eval nproc) --genomeSAindexNbases 13 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${fasta1} ; 
#    cd /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/ ;
done


