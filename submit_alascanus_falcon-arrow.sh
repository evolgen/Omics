#!/bin/bash
#SBATCH --job-name=arrow_alascanus
#
#SBATCH --account=co_genomicdata
#
# QoS: must be savio_long for jobs > 3 days
#SBATCH --qos=savio_lowprio 
#
# Partition:
#SBATCH --partition=savio3_bigmem
#
# Number of tasks needed for use case (example):
#SBATCH -N 1
#SBATCH -n 32
#
# Wall clock limit (7 days in this case):
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON/job_%j.log
#SBATCH --error=/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON/job_%j.err
#
## Command(s) to run (example):

printf "\tRunning alascanus\n"
cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON
mkdir -p arrow && cd arrow
if [ -e "/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq" ]; then
    sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON/alascanus_purged.fasta ~/CONFIGS/falcon/input_bam_alascanus.fofn;
else
    zcat /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq.gz >/clusterfs/genomicdata/rockfish/tmp/alascanus_pacbio.fastq;
    sh ~/SCRIPTS/run_VGP_arrow.sh /clusterfs/genomicdata/rockfish/tmp/alascanus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON/alascanus_purged.fasta ~/CONFIGS/falcon/input_bam_alascanus.fofn;
fi
                  

