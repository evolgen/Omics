#!/usr/bin/bash

set -e

maindir="/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina"

cd ${maindir};
source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh && conda activate mitoz;

for fastq1 in ${maindir}/output/*_*/*-*/corrected/corr.0; do 

    fastq2=$(echo $fastq1 | sed -e 's/\/corr\.0/\/corr.1/');
    workdir=$(echo $fastq1 | sed -e 's/\/corrected\/corr.0$/\/mitogenome/');
    name=$(echo $fastq1 | sed -e 's/\/corrected\/corr.0$//' -e 's/.*\/output\///' -e 's/.*\///')
    mkdir -p ${workdir} &&  cd ${workdir};
    printf "Starting run for ${name} \n"
    if [ ! -e "${workdir}/${name}.result/${name}.fasta" ]; then 
        if [ -e "${workdir}/${name}.tmp/" ]; then
            rm -fr ${workdir}/${name}.tmp/
        fi    
        python3 /global/scratch2/rohitkolora/Software/mitoz/MitoZ.py all2 --run_mode 2 \
            --outprefix ${name} --create_config \
            --fastq1 ${fastq1} --fastq2 ${fastq2} \
            --fastq_read_length 150 --insert_size 450 --fq_size 7 \
            --requiring_taxa 'Chordata' \
            --thread_number 32;
    fi    
    cd ${maindir};

done    


