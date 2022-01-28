#!/usr/bin/bash

filt_vcf=$1

# \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/WTDBG/*/Pacbio/Svim/*/svim.*.filt.vcf >~/RGP/scripts/rearrangements/list_svim_vcffilt.txt2
#\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/Pacbio/Svim/*/svim.*.filt.vcf >~/RGP/scripts/rearrangements/list_svim_vcffilt.txt

module load bcftools

work_dir=$(dirname "$filt_vcf")
reference=$(echo $work_dir | sed -e 's/\/Pacbio\/.*//' -e 's/.*\///')
query=$(echo $work_dir | sed -e 's/.*\/Svim\///' -e 's/\/.*//')

mkdir -p ${work_dir}/BEDPE

printf "  Processing : Reference-$reference\t\tQuery-$query\n";

#bcftools query -i 'SVTYPE=="DEL"' -f '%CHROM\t%POS\t%POS\t%CHROM\t%END\t%END\t%ID\n' $filt_vcf > ${work_dir}/BEDPE/deletions.bedpe
bcftools query -i 'SVTYPE=="DEL"' -f '%CHROM\t%POS\t%END\tSVIM\t%ID\t+\n' $filt_vcf | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/' | awk -v query="$query" -v reference="$reference" 'OFS="\t" {gsub("svim.","",$5); gsub("\\.","_",$5); print $0,".\t.\t.\t.",reference,query}' > ${work_dir}/BEDPE/deletions.bedpe

#bcftools query -i 'SVTYPE=="INV"' -f '%CHROM\t%POS\t%POS\t%CHROM\t%END\t%END\t%ID\n' $filt_vcf > ${work_dir}/BEDPE/inversions.bedpe
bcftools query -i 'SVTYPE=="INV"' -f '%CHROM\t%POS\t%END\tSVIM\t%ID\t+\n' $filt_vcf | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/' | awk -v query="$query" -v reference="$reference" 'OFS="\t" {gsub("svim.","",$5); gsub("\\.","_",$5); print $0,".\t.\t.\t.",reference,query}' > ${work_dir}/BEDPE/inversions.bedpe

#bcftools query -i 'SVTYPE=="INS"' -f '%CHROM\t%POS\t%POS\t%CHROM\t%END\t%END\t%ID\n' $filt_vcf > ${work_dir}/BEDPE/insertions.bedpe
bcftools query -i 'SVTYPE=="INS"' -f '%CHROM\t%POS\t%END\tSVIM\t%ID\t+\n' $filt_vcf | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/' | awk -v query="$query" -v reference="$reference" 'OFS="\t" {gsub("svim.","",$5); gsub("\\.","_",$5); print $0,".\t.\t.\t.",reference,query}' > ${work_dir}/BEDPE/insertions.bedpe

#bcftools query -i 'SVTYPE=="DUP"' -f '%CHROM\t%POS\t%POS\t%CHROM\t%END\t%END\t%ID\n' $filt_vcf > ${work_dir}/BEDPE/duplications.bedpe
bcftools query -i 'SVTYPE=="DUP"' -f '%CHROM\t%POS\t%END\tSVIM\t%ID\t+\n' $filt_vcf | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/' | awk -v query="$query" -v reference="$reference" 'OFS="\t" {gsub("svim.","",$5); gsub("\\.","_",$5); print $0,".\t.\t.\t.",reference,query}' > ${work_dir}/BEDPE/duplications.bedpe

awk -v query="$query" -v reference="$reference" 'OFS="\t" {split($4, f, ";"); split(substr(f[2], 2), g, ":"); print $1, $2, $3, "SVIM", "TRANS_"NR, "+", g[1], g[2], ".\t.",reference,query}' ${work_dir}/candidates/candidates_breakends.bed | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/' > ${work_dir}/BEDPE/translocations.bedpe

paste <(awk 'OFS="\t" {split($4, f, ";"); split(substr(f[2], 2), g, ":"); print $1, $2, $3, "SVIM", "DUP-INT_"NR, "+"}' ${work_dir}/candidates/candidates_int_duplications_source.bed) <(awk -v query="$query" -v reference="$reference" 'OFS="\t" {split($4, f, ";"); split(substr(f[2], 2), g, ":"); print $1, $2, $3, "DUP-INT_"NR,reference,query}' ${work_dir}/candidates/candidates_int_duplications_dest.bed) | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/ '  > ${work_dir}/BEDPE/interspersed_duplications.bedpe

paste <(awk 'OFS="\t" {split($4, f, ";"); split(substr(f[2], 2), g, ":"); print $1, $2, $3, "SVIM", "DUP-TAN_"NR, "+"}' ${work_dir}/candidates/candidates_tan_duplications_source.bed) <(awk -v query="$query" -v reference="$reference" 'OFS="\t" {split($4, f, ";"); split(substr(f[2], 2), g, ":"); print $1, $2, $3, "DUP-TAN_"NR,reference,query}' ${work_dir}/candidates/candidates_tan_duplications_dest.bed) | awk '$1 ~ !/_assoc/ && $7 ~ !/_assoc/'  > ${work_dir}/BEDPE/tandem_duplications.bedpe


