#!/usr/bin/sh

set -e

module load samtools trinity salmon bedtools java jellyfish2 bowtie2

ls /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/star/S_*/*/S*/Aligned.sortedByCoord.out.bam | sed -e 's/\/SRR[0-9]*\/Aligned\.sortedByCoord\.out\.bam//' | sort -u | while read sample; do

  cd ${sample}
  ls ${sample}/S*/Aligned.sortedByCoord.out.bam >${sample}/listed_bams.fofn
  count=$(cat ${sample}/listed_bams.fofn | wc -l)
  
  if [ "$count" -eq 1 ]; then
  unmerged=$(cat ${sample}/listed_bams.fofn | sed -e '/^$/d')
  printf "\n\n#####\tRunning Trinity for ${sample}\t#####\n\n"
  Trinity --seqType fq --max_memory 360G --genome_guided_bam ${unmerged} --genome_guided_max_intron 700000 --CPU 8
  fi

  cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

done

