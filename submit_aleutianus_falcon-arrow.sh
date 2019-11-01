#!/bin/bash
#SBATCH --job-name=arw_aleutianus
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
#SBATCH -n 20
#
# Wall clock limit (7 days in this case):
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON/job_%j.log
#SBATCH --error=/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON/job_%j.err
#
## Command(s) to run (example):

printf "\tRunning aleutianus\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON
mkdir -p arrow && cd arrow
if [ -e "/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/aleutianus_pacbio.fastq" ]; then
    sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/aleutianus_pacbio.fastq /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON/aleutianus_purged.fasta ~/CONFIGS/falcon/input_bam_aleutianus.fofn;
else
    zcat /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/aleutianus_pacbio.fastq.gz >/clusterfs/genomicdata/rockfish/tmp/aleutianus_pacbio.fastq;
    sh ~/SCRIPTS/run_VGP_arrow.sh /clusterfs/genomicdata/rockfish/tmp/aleutianus_pacbio.fastq /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON/aleutianus_purged.fasta ~/CONFIGS/falcon/input_bam_aleutianus.fofn;
fi
                  

