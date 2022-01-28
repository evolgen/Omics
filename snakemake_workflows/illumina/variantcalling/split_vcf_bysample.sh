#!/usr/bin/bash

module load bcftools vcftools
for file in $PWD/../vcfs/seb_aleu_variants.frby.PGA_scaffold_60__9:35000000-35246893.vcf.gz; do
    for sample in `bcftools view -h $file | grep "^#CHROM" | cut -f10-`; do
      vcf_name=$(echo $file | sed -e 's/.*\///' -e 's/\.vcf.gz$//')  
      echo "$sample";
      bcftools view -l1 -Oz -s $sample -o ${vcf_name}.${sample}.vcf.gz $file;
      tabix -p vcf ${vcf_name}.${sample}.vcf.gz $file
    done
done
