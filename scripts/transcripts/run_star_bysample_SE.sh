#!/bin/bash

set -e

module load gcc openmpi trinity bowtie2 jellyfish2 salmon java samtools star

workdir="/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/"

mkdir -p star
cd $workdir

cut -f4 list_samples.txt | sort -u | while read sample; do

  echo "Running for $sample"
  grep -w "$sample" list_samples.txt | sort -k4,4 -k2,2 | while read line1; do
    run=$(echo $line1 | awk -F' ' '{print $1}');
    type=$(echo $line1 | awk -F' ' '{print $2}');
    species=$(echo $line1 | awk -F' ' '{print $3}');
    mkdir -p star/${species}; mkdir -p star/${species}/${sample}; mkdir -p star/${species}/${sample}/${run};

    if [ "$type" != "PAIRED" ]; then
       filesingle=$(ls ./*/${run}.fastq.gz | grep -v adptrm)
       printf "\tStarted singled alignment for $run\n";
       STAR --genomeDir /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/STAR --quantTranscriptomeBan Singleend --runMode alignReads --runThreadN 32 --outSAMtype BAM SortedByCoordinate --readFilesIn $filesingle --outSAMattributes All --outSAMprimaryFlag AllBestScore --outFilterMismatchNmax 10 --readFilesCommand zcat --twopassMode Basic --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFileNamePrefix "./star/${species}/${sample}/${run}"       
    fi

  done

done

