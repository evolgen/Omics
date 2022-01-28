#!/bin/bash

set -e

module load gcc openmpi trinity bowtie2 jellyfish2 salmon java samtools

workdir="/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/"

cd $workdir

for file1 in ${workdir}/Trinity/S_*/*_1.fastq; do

  file2=$(echo $file1 | sed -e 's/_1.fastq$/_2.fastq/')
  file_tmp=$(echo $file1 | sed -e 's/_1.fastq$/_tmp.fastq/')
  filename=$(echo $file1 | sed -e 's/_1.fastq//')
  
  file_F=$(echo $file1 | sed -e 's/_1.fastq$/_Forward.fastq.gz/')
  file_R=$(echo $file1 | sed -e 's/_1.fastq$/_Reverse.fastq.gz/')

  echo "Joining files"
  sed -e 's/^@.* HWI-/@HWI-/' -e '/^@/ s/ length=.*/\/1/' $file1 | gzip -c >${file_F}
  echo "Created $file_F"
  sed -e 's/^@.* HWI-/@HWI-/' -e '/^@/ s/ length=.*/\/2/' $file2 | gzip -c >${file_R}
#  echo "Created $file_R"

  mkdir -p $filename && cd $filename
  Trinity --seqType fq --max_memory 360G --left ${file_F} --right ${file_R} --CPU 8
  cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

done



