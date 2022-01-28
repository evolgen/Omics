#!/usr/bin/bash

filt_vcf=$1

# cat /global/home/users/rohitkolora/RGP/scripts/rearrangements/list_sniffles_vcffilt.txt | parallel -j 1 bash ~/RGP/scripts/rearrangements/sniffles_vcf2bed_byvcf.sh

module load bcftools

work_dir=$(dirname "$filt_vcf")
reference=$(echo $work_dir | sed -e 's/\/Pacbio\/.*//' -e 's/.*\///')
query=$(echo $filt_vcf | sed -e "s/.*sniffles\.${reference}\.//" -e 's/\..*//')
file_name=$(basename "$filt_vcf")

mkdir -p ${work_dir}/BEDPE

printf "  Processing : Reference-$reference\t\tQuery-$query\n";

/global/scratch2/rohitkolora/Software/SURVIVOR/Debug/SURVIVOR vcftobed $filt_vcf 200 -1 ${filt_vcf/.vcf/.bedpe}
awk -F'\t' -v reference="$reference" -v query="$query" 'OFS="\t" {print $1,$2,$5,"SNIFFLES",$11"_"$7,"+",$12,$13,$15,".",reference,query}' ${filt_vcf/.vcf/.bedpe} | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/' >${work_dir}/BEDPE/${file_name/.vcf/.bedpe}

