#!/usr/bin/bash

set -e -o pipefail

coordinates=$@

mkdir -p variants/
module load gcc samtools freebayes java vcftools bcftools

echo "${coordinates}";
if [[ ! -e "variants/seb_aleu_variants.frby.${coordinates}.vcf.gz.tbi" ]]; then
    freebayes --min-mapping-quality 20 --min-base-quality 20 \
    -r ${coordinates} \
    -f /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/fSebaleuF.fasta \
    --min-coverage 5 --gvcf -g 1000 \
    -L rg_bamfile_list.txt |
    vcffilter -f "QUAL > 20" >variants/seb_aleu_variants.frby.${coordinates}.vcf;
bgzip variants/seb_aleu_variants.frby.${coordinates}.vcf;
tabix -f -p vcf variants/seb_aleu_variants.frby.${coordinates}.vcf.gz &

fi

