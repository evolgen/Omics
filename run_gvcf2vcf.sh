#!/usr/bin/bash

set -e

module load vcftools samtools bcftools freebayes gcc;

file1=$@;
file2=$(echo $file1 | sed -e 's/.vcf.gz$/.VCF.gz/');

echo "Converting $file1"; 
if [ ! -e "${file2}.tbi" ]; then
    bcftools convert --gvcf2vcf $file1 \
        -f /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/WTDBG2/arrow/frby_polish/salsa_dnase2/scaffolds/scaffolds_FINAL.fasta | bcftools view --exclude-types indels | vcf-sort | \
        bgzip -c >$file2;
    tabix -p vcf $file2;
fi


