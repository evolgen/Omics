#!/usr/bin/bash

set -e

#####
###sh ~/SCRIPTS/run_raconpolish_illumina.sh \
###		/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/Seb10_S144.fastq \
###		/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/wtdbg2_007.V1.ctg.fa
#####

if [ $# -ne 2 ]; then
	printf "\n\n\tPlease provide 3 parameters as follows\n";
	printf "\t\trun_raconpolishillumina.sh FASTQ REF_FASTA\n\n";
	exit 1;
fi

fastq=$1
reference=$2

module load racon minimap2 samtools

count=0
ln -sf $reference ${reference}.${count}.fasta
samtools faidx ${reference}.${count}.fasta

for iter in `seq 1 4`; do
	version=$((count+1));
	if [ ! -e "map_${reference}.${count}.sam" ]; then
		minimap2 -t 32 -ax sr ${reference}.${count}.fasta ${fastq} >map_${reference}.${count}.sam;
	fi
	racon -u -e 0.2 -t 32 ${fastq} map_${reference}.${count}.sam ${reference}.${count}.fasta >${reference}.${version}.fasta;
	samtools faidx ${reference}.${version}.fasta;
	rm map_${reference}.${count}.sam;
	count=$((count+1));
done



