#!/bin/bash
#SBATCH --job-name=rnaseq_SE
#
#SBATCH --account=fc_genomicdata
#
# QoS: must be savio_long for jobs > 3 days
#
# Partition:
#SBATCH --partition=savio
#
# Number of tasks needed for use case (example):
#SBATCH --ntasks=20
#SBATCH --nodes=1
#
# Wall clock limit (7 days in this case):
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/job_rnaseq_SE.log
#SBATCH --error=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/job_rnaseq_SE.err
#
## Command(s) to run (example):

set -e

module load gcc openmpi fastp flash

cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

for file in ./S*/*[!_][0-9].fastq.gz; do

  fileout=$(basename "$file")
  filename=$(echo $file | sed -e 's/.fastq.gz//')

  fastp --trim_poly_g --overrepresentation_analysis --html ${filename}.html -j ${filename}.json --thread 20 --length_required 35 --overlap_diff_limit 2 --trim_poly_x --disable_quality_filtering --n_base_limit 1 --in1 $file --out1 ${file}.adptrm

done

