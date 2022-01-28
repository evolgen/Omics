#!/usr/bin/bash

filt_in=$1

# cat ~/RGP/scripts/rearrangements/list_sedef_bedfilt.txt | parallel -j 1 bash ~/RGP/scripts/rearrangements/sedef_bed_convert.sh

module load bcftools

work_dir=$(dirname "$filt_in")
reference=$(echo $work_dir | sed -e 's/\/SEDEF.*//' -e 's/.*\///')
file_name=$(basename "$filt_in")

mkdir -p ${work_dir}/BEDPE

printf "  Processing : Reference-$reference\t\tQuery-$reference\n";

awk -F'\t' -v reference="$reference" 'OFS="\t" {if($1~!/^Sebast/) {print reference"."$1, $2, $3, "SEDEF", "SEG-DUP_"NR, $9, reference"."$4, $5,$6,$NF, reference, reference} else {print $1, $2, $3, "SEDEF", "SEG-DUP_"NR, $9, $4, $5,$6,$NF, reference, reference} }' $filt_in >${work_dir}/BEDPE/sedef.bedpe

