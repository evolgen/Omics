#!/usr/bin/bash

set -e

#####
###sh ~/SCRIPTS/run_quiverpolish.sh \
###		/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/umbrosus_pacbio.fastq \
###		/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/wtdbg2_007.V1.ctg.fa
#####

if [ $# -ne 2 ]; then
	printf "\n\n\tPlease provide 2 parameters as follows\n";
	printf "\t\trun_quiverpolish.sh FASTQ REF_FASTA\n\n";
	exit 1;
fi

xml=$1
reference=$2

source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh
conda activate py2.7 && conda activate pb-assembly
module load minimap2 samtools #smrtlink

count=0
ln -sf $reference ${reference}.${count}.fasta
samtools faidx ${reference}.${count}.fasta

for iter in `seq 1 5`; do
	version=$((count+1));
	if [ ! -e "pbmap_${reference}.${count}.bam" ]; then
		pbmm2 align --preset SUBREAD ${reference}.${count}.fasta $xml pbmap_${reference}.${count}.bam --sort -j 26 -J 6;
	fi	
	if [ ! -e "pbmap_${reference}.${count}.bam.pbi" ]; then
		pbindex pbmap_${reference}.${count}.bam;
	fi	
	variantCaller pbmap_${reference}.${count}.bam --algorithm=arrow -j 32 -r ${reference}.${count}.fasta -o ${reference}.${version}.fasta;
	samtools faidx ${reference}.${version}.fasta;
	count=$((count+1));
done



