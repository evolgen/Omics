#!/usr/bin/bash

set -e

module load wtdbg2 minimap2 samtools/1.8

wtdbg2 -p 21 -K 800.010010 -A -S 1.000000 -s 0.050000 -g 1g -e 3 -L 5000 -X 50 -o ruberrimus_wtdbg2_001.V0 -i /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/ruberrimus_pacbio.fastq.gz 1>ruberrimus_wtdbg2_001.LOG 2>>ruberrimus_wtdbg2_001.LOG

~/local_modules_sw/wtdbg2/2.3/wtpoa-cns -t 32 -i ruberrimus_wtdbg2_001.V0.ctg.lay.gz -fo ruberrimus_wtdbg2_001.V0.ctg.fa 1>>ruberrimus_wtdbg2_001.LOG 2>>ruberrimus_wtdbg2_001.LOG
gzip < ruberrimus_wtdbg2_001.V0.ctg.fa >ruberrimus_wtdbg2_001.V0.ctg.fa.gz

for count in `seq 0 5`; do

  version=$((count+1))
  ~/local_modules_sw/minimap2/2.15/minimap2 -t 32 -k 19 -w 10 -a ruberrimus_wtdbg2_001.V${count}.ctg.fa /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/ruberrimus_canu.correctedReads.fasta.gz >sam.sam;
  samtools view -Sb sam.sam >bam.bam; 
  rm sam.sam;
  samtools sort -o ruberrimus_wtdbg2_001.V${count}.ctg.map.srt.bam bam.bam;
  rm bam.bam;
  ~/local_modules_sw/wtdbg2/2.3/wtpoa-cns -t 32 -d ruberrimus_wtdbg2_001.V${count}.ctg.fa -i ruberrimus_wtdbg2_001.V${count}.ctg.map.srt.bam -fo ruberrimus_wtdbg2_001.V${version}.ctg.fa 1>>ruberrimus_wtdbg2_001.LOG 2>>ruberrimus_wtdbg2_001.LOG; 
  gzip < ruberrimus_wtdbg2_001.V${version}.ctg.fa >ruberrimus_wtdbg2_001.V${version}.ctg.fa.gz;

done


