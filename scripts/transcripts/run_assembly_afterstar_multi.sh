#!/usr/bin/sh

set -e

module load samtools trinity salmon bedtools java jellyfish2 bowtie2

ls /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/star/S_*/*/*/Aligned.sortedByCoord.out.bam | sed -e 's/\/SRR[0-9]*\/Aligned\.sortedByCoord\.out\.bam//' | sort -u | while read sample; do

  cd ${sample}
  ls ${sample}/S*/Aligned.sortedByCoord.out.bam >${sample}/listed_bams.fofn
  count=$(cat ${sample}/listed_bams.fofn | wc -l)

  if [ "$count" -gt 1 ]; then
  printf "Merging \t ${sample}\n"
  files_list=$(cat ${sample}/listed_bams.fofn | tr "\n" " " | sed -e 's/$/\n/')

#  samtools merge -@ 32 -fcp ${sample}/merged.bam ${files_list}
  printf "\n\n#####\tRunning Trinity for ${sample}\t#####\n\n"
  Trinity --seqType fq --max_memory 360G --genome_guided_bam ${sample}/merged.bam --genome_guided_max_intron 700000 --CPU 8
  fi

  cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

done

