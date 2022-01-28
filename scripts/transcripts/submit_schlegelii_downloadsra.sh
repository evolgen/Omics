#!/bin/bash
#SBATCH --job-name=schlegelii
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
#SBATCH --output=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/S_schlegelii/job_%j.log
#SBATCH --error=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/S_schlegelii/job_%j.err
#
## Command(s) to run (example):

module load gcc openmpi sratoolkit/2.9.4

cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/S_schlegelii
fasterq-dump SRR4409372 --split-3 --skip-technical -e 20 -m 1000MB
fasterq-dump SRR4409389 --split-3 --skip-technical -e 20 -m 1000MB
fasterq-dump SRR4409390 --split-3 --skip-technical -e 20 -m 1000MB
