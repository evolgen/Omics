#!/usr/bin/sh

set -e

name=$1
seq=$@

pacbioseq=$(echo $seq | sed -e 's/*.\.sh//' -e "s/.*${name} //" -e 's/^/-i /' -e 's/ /-i /g')

module load wtdbg2 minimap2 samtools/1.8

 wtdbg2 -t 32 -x corrected -S 2 -o $name $pacbioseq 1>$name.LOG 2>>$name.ERR
 
 ~/local_modules_sw/wtdbg2/2.3/wtpoa-cns -t 32 -i $name.ctg.lay.gz -fo $name.V0.ctg.fa 1>>$name.LOG 2>>$name.ERR; 
 gzip < $name.V0.ctg.fa >$name.V0.ctg.fa.gz

for count in `seq 0 5`; do

  version=$((count+1))
  ~/local_modules_sw/minimap2/2.15/minimap2 -t 32 -x map-pb -a $name.V${count}.ctg.fa /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/SEB10_Cell1/m54050R1_181211_022522.subreads.fastq.gz /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/SEB10_Cell2/m54050R1_181212_232858.subreads.fastq.gz /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/SEB10_Cell3/m54050R1_181213_094204.subreads.fastq.gz /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/SEB10_Cell4/m54290_190124_022926.subreads.fastq.gz | samtools view -Sb - >$name.ctg.map.bam;
  samtools sort -@ 32 -o $name.V${count}.ctg.map.srt.bam $name.ctg.map.bam; 
  rm $name.ctg.map.bam;
  samtools view $name.V${count}.ctg.map.srt.bam | ~/local_modules_sw/wtdbg2/2.3/wtpoa-cns -t 32 -d $name.V${count}.ctg.fa -i - -fo $name.V${version}.ctg.fa 1>>$name.LOG 2>>$name.ERR;
  gzip < $name.V${version}.ctg.fa >$name.V${version}.ctg.fa.gz;

done




