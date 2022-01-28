#!/bin/bash
#SBATCH --job-name=rnaseq_PE
#
#SBATCH --account=fc_genomicdata
#
# QoS: must be savio_long for jobs > 3 days
#
# Partition:
#SBATCH --partition=savio
#
# Number of tasks needed for use case (example):
#SBATCH --ntasks=32
#SBATCH --nodes=1
#
# Wall clock limit (7 days in this case):
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/job_rnaseq_PE.log
#SBATCH --error=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/job_rnaseq_PE.err
#
## Command(s) to run (example):

set -e

module load gcc openmpi trinity bowtie2 jellyfish2 salmon java samtools

cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

for file1 in ./Trinity/S_*/*_1.fastq; do

  file2=$(echo $file1 | sed -e 's/_1.fastq$/_2.fastq/')
  file_tmp=$(echo $file1 | sed -e 's/_1.fastq$/_tmp.fastq/')
  filename=$(echo $file1 | sed -e 's/_1.fastq//')
  
  file_F=$(echo $file1 | sed -e 's/_1.fastq$/_Forward.fastq.gz/')
  file_R=$(echo $file1 | sed -e 's/_1.fastq$/_Reverse.fastq.gz/')

  cat $file1 $file_tmp | gzip -c >$file_F
  cat $file2 | gzip -c >$file_R

  mkdir $filename && cd $filename
  Trinity --seqType fq --max_memory 360G --left $file_F --right $file_R --CPU 8
  cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

done



