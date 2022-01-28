#!/usr/bin/bash

set -e -o pipefail

coordinates=$@

mkdir -p vcfs/
module load gcc samtools freebayes java vcftools

echo "${coordinates}";
if [[ ! -e "vcfs/seb_aleu_variants.frby.${coordinates}.vcf.gz.tbi" ]]; then
    freebayes --min-mapping-quality 20 --min-base-quality 20 \
    -r ${coordinates} \
    -f /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/fSebaleuF.fasta \
    --report-monomorphic --min-coverage 5 \
    -L rg_bamfile_list.txt >vcfs/seb_aleu_variants.frby.${coordinates}.vcf;
bgzip vcfs/seb_aleu_variants.frby.${coordinates}.vcf;
tabix -f -p vcf vcfs/seb_aleu_variants.frby.${coordinates}.vcf.gz &

fi

