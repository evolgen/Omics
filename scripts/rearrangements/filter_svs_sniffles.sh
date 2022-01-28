#!/usr/bin/bash

set -e

#### Filter SVIM vcf on a quality of ten

vcf=$1
name=$2

if [ "$vcf" ~ ".vcf.gz$" ]; then
    zcat $vcf | 
        grep -v -e 'IMPRECISE;' |
        sed -e "/^#CHROM/ s/\/global\/.*/${name}/" >${vcf/%.vcf.gz/.filt.vcf}
else        
    grep -v -e 'IMPRECISE;' $vcf |
        sed -e "/^#CHROM/ s/\/global\/.*/${name}/" >${vcf/%.vcf/.filt.vcf}
fi        


