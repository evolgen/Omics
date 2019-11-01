#!/bin/bash
#SBATCH --job-name=CANU_alascanus
#
#SBATCH --account=co_genomicdata
#
# QoS: must be savio_long for jobs > 3 days
#SBATCH --qos=savio_lowprio 
#
# Partition:
#SBATCH --partition=savio3_xlmem
#
# Number of tasks needed for use case (example):
#SBATCH -N 1
#SBATCH -n 32
#
# Wall clock limit (7 days in this case):
#SBATCH --time=3-00:00:00
#SBATCH --output=/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_alascanus/05.28.2019/Assembly/canu/job_%j.log
#SBATCH --error=/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_alascanus/05.28.2019/Assembly/canu/job_%j.err
#
## Command(s) to run (example):

source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh
conda activate py2.7
module load canu samtools minimap2 wtdbg2 java

cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_alascanus/05.28.2019/Assembly/canu

canu merylConcurrency=10 hapConcurrency=10 cormhapConcurrency=10 obtovlConcurrency=10 utgovlConcurrency=10 corConcurrency=10 ovbConcurrency=10 ovsConcurrency=10 redConcurrency=10 oeaConcurrency=10 batConcurrency=10 cnsConcurrency=10 corMinCoverage=0 batThreads=64 maxThreads=64 batMemory=1550g minMemory=1500g maxMemory=1560g corOutCoverage=150 correctedErrorRate=0.105 OvlMerThreshold=500 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" -p alascanus_canu mhapSensitivity=normal -d run1 rawErrorRate=0.36 genomeSize=0.95g gridOptions="--time=72:00:00 --account=co_genomicdata --qos=savio_lowprio --partition=savio3_bigmem -N 1 -n 32" useGrid=false -pacbio-raw /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq.gz 1>>LOG2 2>>LOG2


