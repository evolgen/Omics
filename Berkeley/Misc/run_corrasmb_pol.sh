#!/usr/bin/sh

set -e

name=$1
###seq=$@
corrpacbioseq=$2
pacbioseq=$3

###pacbioseq=$(echo $seq | sed -e 's/*.\.sh//' -e "s/.*${name} //" -e 's/^/-i /' -e 's/ /-i /g')

module load wtdbg2 minimap2 samtools/1.8

printf "\n\tAssembling with WTDBG2 - ${name} \n";

if [[ ! -e "${name}.V0.ctg.fa.gz" ]]; then
     wtdbg2 -t 32 -x corrected -g 1g -e 3 -L 5000 -X 45 -S 2 -o ${name} ${corrpacbioseq} 1>${name}.LOG 2>${name}.ERR
     ~/local_modules_sw/wtdbg2/2.3/wtpoa-cns -t 32 -i ${name}.ctg.lay.gz -fo ${name}.V0.ctg.fa 1>>${name}.LOG 2>>${name}.ERR; 
     gzip < ${name}.V0.ctg.fa >${name}.V0.ctg.fa.gz
fi
     

for count in `seq 0 5`; do

  version=$((count+1))
  ~/local_modules_sw/minimap2/2.15/minimap2 -t 32 -k 19 -w 10 -a ${name}.V${count}.ctg.fa ${pacbioseq} | \
      samtools view -bS - >${name}.ctg.map.bam;
  samtools sort -T ./ -m 1G -@ 32 -o ${name}.V${count}.ctg.map.srt.bam ${name}.ctg.map.bam; 
  rm ${name}.ctg.map.bam;
  samtools view ${name}.V${count}.ctg.map.srt.bam | \
      ~/local_modules_sw/wtdbg2/2.3/wtpoa-cns -t 32 -d ${name}.V${count}.ctg.fa -i - -fo ${name}.V${version}.ctg.fa \
      1>>${name}.LOG 2>>${name}.ERR;
  gzip < ${name}.V${version}.ctg.fa >${name}.V${version}.ctg.fa.gz;

done




