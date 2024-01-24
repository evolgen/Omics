#!/usr/bin/bash

set -e

module load vcftools samtools bcftools freebayes gcc;

file1=$@;
file2=$(echo $file1 | sed -e 's/Dedup.bam$/Variant.snp.vcf.gz/');
sample=$(echo $file1 | sed -e 's/\/calls//' -e 's/.*\/S-/S-/' -e 's/.*\/B-/B-/' -e 's/.*\/H-/H-/' -e 's/\/.*//')

echo "Converting $file1 using $sample"; 
if [ ! -e "${file2}.tbi" ]; then
    freebayes -f /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/WTDBG2/arrow/frby_polish/salsa_dnase2/scaffolds/scaffolds_FINAL.fasta --min-coverage 8 --exclude-unobserved-genotypes --use-mapping-quality -m 20 -q 10 -i -X -u ${file1} |
        vcffilter -f "QUAL > 20" |
        sed -e "0,/^#CHROM/ s/unknown$/${sample}/" |
        bgzip -c >${file2};
    tabix -p vcf $file2;
fi


