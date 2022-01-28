#!/usr/bin/bash

set -e

#### Filter SVIM vcf on a quality of ten

vcf=$1

if [[ "$vcf" =~ ".vcf.gz$" ]]; then
    zcat $vcf | 
        awk '{ if($1 ~ /^#/) { print $0 } else { if($6>=10) { print $0 } } }' >${vcf/%.vcf.gz/.filt.vcf}
else        
    awk '{ if($1 ~ /^#/) { print $0 } else { if($6>=10) { print $0 } } }' $vcf >${vcf/%.vcf/.filt.vcf}
fi        


