#!/usr/bin/bash

set -e

##### Polishing reference with Illumina data using Freebayes #####
# author - Rohit Kolora #
# date - 14th August, 2019 #
#####
### sh ~/freebayes_polish.sh \
###		/global/scratch2/rohitkolora/Illumina_R1.fastq.gz \
###		/global/scratch2/rohitkolora/Illumina_R2.fastq.gz \
###		/global/scratch2/rohitkolora/reference.fasta
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

module load minimap2 samtools java freebayes

bam=$(echo $reference | sed -e 's/.*\///' -e 's/\.fasta/.1.bam/' -e 's/\.fas/.1.bam/' -e 's/\.fa/.1.bam/')
bam2=$(echo $reference | sed -e 's/.*\///' -e 's/\.fasta/.2.bam/' -e 's/\.fas/.2.bam/' -e 's/\.fa/.2.bam/')
consensus=$(echo $reference | sed -e 's/.*\///' -e 's/\.fa$/.frby/' -e 's/\.fas$/.frby/' -e 's/\.fasta$/.frby/')
final_consen=$(echo ${consensus}.2.fasta | sed -e 's/.t1.frby.2.fasta/.t2.fasta/')

if [ ! -e "${consensus}.2.fasta.fai" ]; then     ### MAIN

if [ ! -e "${fasta}.fai" ]; then                ### INDEXING
    samtools faidx ${fasta}
fi    

printf "\tMapping short reads for 1st cycle\n";     ### 1st Mapping
if [ ! -e "${bam}.bai" ]; then
    if [ ! -e "${reference}.sort.bam.bai" ]; then
        minimap2 -t 24 -ax sr ${fasta} ${fastq1} ${fastq2} | samtools view -bh >${reference}.bam;
        printf "\tSorting mapped reads\n";
        samtools sort -@ 24 -T ./ -o ${reference}.sort.bam ${reference}.bam; rm ${reference}.bam &
    fi

    printf "\tRunning GATK for duplicates for 1st round\n";
    gatk --java-options "-Xmx100G" MarkDuplicates --TMP_DIR $PWD --INPUT ${reference}.sort.bam --OUTPUT ${bam} --METRICS_FILE ${bam}.metrics --MAX_FILE_HANDLES 15000;
    samtools index ${bam}
    rm ${reference}.sort.bam &
    conda deactivate 
fi

module load samtools freebayes vcftools gcc bcftools

printf "\tFreebayesing now for 1st iteration\n";    ### 1st Freebayesing
if [ ! -e "${bam/%.bam/.vcf.gz}.tbi" ]; then
    freebayes --limit-coverage 200 -f ${fasta} ${bam} | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f ${fasta} -o ${bam/%.bam/.bcf} --threads 32;
    bcftools index ${bam/%.bam/.bcf};
    bcftools consensus -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $fasta ${bam/%.bam/.bcf} > ${consensus}.1.fasta 
    bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=32 ${bam/%.bam/.bcf} > ${bam/%.bam/.vcf.gz}
    tabix -p vcf -f ${bam/%.bam/.vcf.gz}
fi

if [ ! -e "${bam2}.bai" ]; then
    if [ ! -e "${reference}.sort2.bam.bai" ]; then
        printf "\t\t\nROUND 2\n\n\tMapping short reads\n";  ### 2nd Mapping
        minimap2 -t 24 -ax sr ${consensus}.1.fasta ${fastq1} ${fastq2} | samtools view -bh >${reference}.bam;
        printf "\tSorting mapped reads\n";
        samtools sort -@ 24 -T ./ -o ${reference}.sort2.bam ${reference}.bam; rm ${reference}.bam &
    fi

    
    printf "\tRunning GATK for duplicates - round 2\n";
    gatk --java-options "-Xmx100G" MarkDuplicates --TMP_DIR $PWD --INPUT ${reference}.sort2.bam --OUTPUT ${bam2} --METRICS_FILE ${bam2}.metrics --MAX_FILE_HANDLES 15000;
    samtools index ${bam2} && samtools faidx ${consensus}.1.fasta;
    rm ${reference}.sort2.bam &
    conda deactivate
fi

if [ ! -e "${bam2/%.bam/.vcf.gz}.tbi" ]; then
    printf "\tFreebayesing now - round 2\n";            ### 2nd Freebayesing
    freebayes --limit-coverage 200 -f ${consensus}.1.fasta ${bam2} | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f ${consensus}.1.fasta -o ${bam2/%.bam/.bcf} --threads 24;
    bcftools index ${bam2/%.bam/.bcf};
    bcftools consensus -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f ${consensus}.1.fasta ${bam2/%.bam/.bcf} > ${consensus}.2.fasta
    bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=32 ${bam2/%.bam/.bcf} > ${bam2/%.bam/.vcf.gz}
    tabix -p vcf -f ${bam2/%.bam/.vcf.gz}
fi

ln -sf ${consensus}.2.fasta ${final_consen}         ### FINAL consensus

fi          ###MAIN


