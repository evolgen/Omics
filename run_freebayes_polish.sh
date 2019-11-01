#!/usr/bin/bash

set -e

#####
###sh ~/SCRIPTS/run_freebayes_polish.sh \
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
fastq2=$2
fasta=$3

reference=$(echo $fasta | sed -e 's/.*\///')

module load minimap2 samtools java

bam=$(echo $reference | sed -e 's/.*\///' -e 's/\.fasta/.bam/' -e 's/\.fas/.bam/' -e 's/\.fa/.bam/')
consensus=$(echo $reference | sed -e 's/.*\///' -e 's/\.fa$/.frby.fasta/' -e 's/\.fas$/.frby.fasta/' -e 's/\.fasta$/.frby.fasta/')

printf "\tMapping short reads\n";
minimap2 -t 24 -ax sr ${fasta} ${fastq1} ${fastq2} | samtools view -bh >${reference}.bam;
printf "\tSorting mapped reads\n";
samtools sort -@ 24 -T ./ -o ${reference}.sort.bam ${reference}.bam; rm ${reference}.bam &

source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh
conda activate gatk

printf "\tRunning GATK for duplicates\n";
/global/scratch2/rohitkolora/Software/gatk-4.1.2/gatk --java-options "-Xmx100G" MarkDuplicates --TMP_DIR $PWD --INPUT ${reference}.sort.bam --OUTPUT ${bam} --METRICS_FILE ${bam}.metrics --MAX_FILE_HANDLES 15000;
rm ${reference}.sort.bam &
conda deactivate 

module load samtools freebayes vcftools gcc
samtools index ${bam};
printf "\tFreebayesing now\n";
freebayes -X -u --min-mapping-quality 30 --min-alternate-count 8 --min-coverage 8 --max-coverage 500 -f ${fasta} ${bam} | vcffilter -f "QUAL > 39 & AF > 0.9" | vcffilter -f "TYPE = snp | TYPE = ins | TYPE = del" >${bam/%.bam/.vcf};

printf "\tCreating index of vcf now\n";
bgzip ${bam/%.bam/.vcf}
tabix -p vcf ${bam/%.bam/.vcf.gz}

printf "\tConsensusing at the moment\n";
cat ${fasta} | vcf-consensus ${bam/%.bam/.vcf.gz} > ${consensus}
printf "All Done consensus making\n";

