#!/usr/bin/bash

set -e

if [ $# -ne 3 ]; then
    printf "\t Please run as follows -\n";
    printf "\t\t sh ~/RGP/run_polishing_freebayes.sh /PATH/TO/REF/FASTA /PATH/TO/PAIREND/FASTQ1 /PATH/TO/PAIREND/FASTQ2\n\n";
    exit 1;
fi    

reference_genome=$1
pe_fq1=$2
pe_fq2=$3

directory1=$(dirname "${reference_genome}")
sample=$(echo ${reference_genome} | sed -e 's/.*\///' -e 's/\..*//' )
workdir=$(echo $directory1 | sed -e 's/$/\/frby_polish/')
mkdir -p $workdir && cd $workdir

CPU=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
CPU2=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l | awk '{print $0/2}')
MEM=$(free | awk 'NR==2 {print $4/1000000}' | sed -e 's/\.[0-9]*//')

source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh;
module load minimap2 samtools bcftools freebayes gcc;

ln -sf ${reference_genome} . 
#if [[ ! -e "${reference_genome}.fai" ]]; then
samtools faidx ${reference_genome} &
#fi


fasta="${reference_genome}";
bam="aligned.bam";
bcf1="output1.bcf";

for count in 2 3; do

    printf "ROUND $count\n";
    cou=$((count-1));
#    if [[ "$count" -eq 3 ]]; then
    printf "\tRunning Minimap2 and sorting bam simultaneously\n";
#    minimap2 -t ${CPU2} -ax sr ${reference_genome} \
    minimap2 -t ${CPU2} -ax sr ${sample}.t${cou}.fasta \
        ${pe_fq1} ${pe_fq2} | \
        samtools view -h - | \
        samtools sort -@ ${CPU2} -m 2G \
        -T /clusterfs/genomicdata/rockfish/tmp \
        -o output.bam;
    printf "\tIndexing the mapped bam and running GATK markduplicates\n";
    samtools index output.bam;
#    fi
    conda activate gatk;
    /global/scratch2/rohitkolora/Software/gatk-4.1.2/gatk \
            --java-options "-Xmx${MEM}G" MarkDuplicates \
            --TMP_DIR /clusterfs/genomicdata/rockfish/tmp \
            --INPUT output.bam \
            --OUTPUT aligned.${count}.bam \
            --METRICS_FILE markdup.metrics \
            --MAX_FILE_HANDLES 15000;
    conda deactivate;
    printf "\tIndexing the processed bam and removing previous bam\n";
    rm output.bam output.bam.bai &
    samtools index aligned.${count}.bam;

    printf "\tRunning Freebayes \n";
    freebayes --bam aligned.${count}.bam --max-coverage 200 -f ${sample}.t${cou}.fasta | \
        bcftools view --no-version -Ou > ${bcf1};

    printf "\tRunning bcftools normalisation \n";
    bcftools view -Ou -e'type="ref"' ${bcf1} | \
        bcftools norm -Ob -f ${sample}.t${cou}.fasta \
        -o ${sample}.${count}.bcf --threads $((CPU-1)) \
        && rm ${bcf1};
    printf "\tRunning bcftools to generate variant file \n";
    bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz \
        --threads=${CPU} ${sample}.${count}.bcf \
        >${sample}.${count}.changes.vcf.gz;
    printf "\tIndexing the variant file for post processing\n";
    bcftools index ${sample}.${count}.changes.vcf.gz; 

    printf "\tGenerating the consensus\n";
    bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla \
        -f ${sample}.t${cou}.fasta ${sample}.${count}.changes.vcf.gz \
        > ${sample}.t${count}.fasta;
    samtools faidx ${sample}.t${count}.fasta;
    printf "\tCounting the number of changes for QV\n";
    bcftools view -H -i 'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ov ${sample}.${count}.changes.vcf.gz | \
        awk -F "\t" '{print $4"\t"$5}' | \
        awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' \
            > ${sample}.${count}.numvar
    printf "\nRound $count\t:\tNum. bases affected: `cat ${sample}.${count}.numvar`\n"

done



