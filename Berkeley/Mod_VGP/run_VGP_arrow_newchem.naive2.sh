#!/usr/bin/bash

set -e

#####
###sh ~/SCRIPTS/run_VGP_arrow.sh \
###		/path/to/input_pacbio.fastq \
###		/path/to/ref.fa \
###     /path/to/input_bam.fofn
#####

if [ $# -ne 3 ]; then
	printf "\n\n\tPlease provide 3 parameters as follows\n";
	printf "\t\trun_VGP_arrow_newchem.naive.sh FASTQ REF_FASTA BAM_FOFN\n\n";
	exit 1;
fi

fastq=$1
fastaseq=$2
bam_fofn=$3

reference=$(echo $fastaseq | sed -e 's/.*\///' -e 's/\.fasta$//' -e 's/\.fas$//' -e 's/\.fa$//' -e 's/\.FASTA$//' -e 's/\.FAS$//' -e 's/\.FA$//')

#source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh
#conda activate pb-assembly
module load minimap2 samtools #smrtlink

count=0
mkdir -p arrow && cd arrow

printf "\tCreating reference file\n";
ln -sf $fastaseq ${reference}.${count}.fasta
samtools faidx ${reference}.${count}.fasta

uniqID=$(echo $BASHPID)

for iter in 1; do
#    conda activate smrtlink
    pbindex pbmap_${reference}.pbsort.bam ;
    echo "  Starting consensus caller" ;
    variantCaller pbmap_${reference}.pbsort.bam --algorithm=arrow -j 32 -r ${reference}.${count}.fasta -o ${reference}.t1.fasta;
	samtools faidx ${reference}.t1.fasta;
	count=$((count+1));
done

printf "\tPolished to produce ${reference}.t1.fasta\n";

