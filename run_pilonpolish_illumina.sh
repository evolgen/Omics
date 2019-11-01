#!/usr/bin/bash

set -e

#####
###sh ~/SCRIPTS/run_pilonpolish_illumina.sh \
###		/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_entomelas/S-entomelas_SEB-8/Seb8_S126_R1_001.fastq.gz \
###		/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_entomelas/S-entomelas_SEB-8/Seb8_S126_R2_001.fastq.gz \
###		/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/polish/wtdbg2_002.V2.ctg.fa
#####

if [ $# -ne 3 ]; then
	printf "\n\n\tPlease provide 3 parameters as follows\n";
	printf "\t\trun_raconpolishillumina.sh FASTQ1 FASTQ2 REF_FASTA\n\n";
	exit 1;
fi

fastq1=$1
fastq1=$2
ref_seq=$3
reference=$(echo $ref_seq | sed -e 's/.*\///' -e 's/\.fasta$//' -e 's/\.fas$//' -e 's/\.fa$//')

module load java minimap2 samtools pilon

count=0
ln -sf $ref_seq ${reference}.pilon.${count}.fasta
printf "\tIndexing the reference\n";
samtools faidx ${reference}.pilon.${count}.fasta

for iter in `seq 1 5`; do

	version=$((count+1));
	bamname=$(echo "${reference}" | sed -e 's/\.fasta$//' -e 's/\.fas$//' -e 's/\.fa$//' -e 's/.*\///')
	
	if [ ! -e "$sort_${bamname}.bam" ]; then
		printf "\tMapping illumina reads with Minimap2\n";
		minimap2 -t 32 -ax sr ${reference}.pilon.${count}.fasta ${fastq1} ${fastq2} | samtools view -bSh - >${bamname}.${count}.bam;

		printf "\tSorting and indexing mapped file\n";
		samtools sort -@ 32 -T $PWD -o sort_${bamname}.${count}.bam ${bamname}.${count}.bam; 
		samtools index sort_${bamname}.${count}.bam;
		rm ${bamname}.${count}.bam &
	fi

	printf "\tExtracting only mapped reads\n";
	samtools view -hb -f 2 sort_${bamname}.${count}.bam >paired.bam;
	samtools view -hb -F 8 -f 4 sort_${bamname}.${count}.bam >matemap.bam;
	samtools view -hb -f 8 -F 4 sort_${bamname}.${count}.bam >readmap.bam;
	
	printf "\tMerging unpaired map files\n";
	samtools merge both.bam readmap.bam matemap.bam;
	samtools sort -@ 32 -T $PWD -o single.bam both.bam;
	rm both.bam readmap.bam matemap.bam &
	printf "\tIndexing mapped single bams\n";
	samtools index paired.bam; samtools index single.bam;
	
	printf "\tRunning PILON\n";
	java -Xmx1550g -jar ~/local_modules_sw/pilon/pilon.jar --genome ${reference}.pilon.${count}.fasta --frags paired.bam --unpaired single.bam --output ${reference}.pilon.${version} --outdir $PWD --diploid --chunksize 200000000 1>>LOG 2>>LOG;
	
	rm paired.bam single.bam & 
	printf "\tIndexing the pilon output fasta\n";
	samtools faidx ${reference}.pilon.${version}.fasta;
#	rm sort_${bamname}.bam;
	count=$((count+1));
	printf "\n\t\tALL DONE :D\n\n"
done



