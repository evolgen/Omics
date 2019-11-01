#!/usr/bin/bash

set -e

module load vcftools samtools bcftools freebayes gcc;

file1=$@;
file1_2=$(echo $file1 | sed -e 's/.vcf.gz$/.tmp.vcf.gz/');
file2=$(echo $file1 | sed -e 's/.vcf.gz$/.SNPs.vcf.gz/');

echo "Converting $file1"; 
if [ ! -e "${file2}.tbi" ]; then
    vcftools --gzvcf $file1 --remove-indels --min-alleles 2 \
        --recode --recode-INFO-all \
        --out SNPs_only --stdout | bgzip -c > ${file1_2};
    vcf-sort ${file1_2} | bgzip -c >$file2 && rm ${file1_2};
    tabix -p vcf $file2;
fi


