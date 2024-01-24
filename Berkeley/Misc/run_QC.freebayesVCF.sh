#!/usr/bin/bash

set -e

if [ $# -ne 1 ]; then
    printf "\t Please run as follows -\n";
    printf "\t\t sh ~/RGP/run_QC.freebayesVCF.sh /PATH/TO/ARROW/FRBY_POLISH/VCF.GZ\n\n";
    exit 1;
fi    

reference_vcf=$1

directory1=$(dirname "${reference_vcf}")
sample=$(echo ${reference_vcf} | sed -e 's/.*\///' -e 's/\..*//' )
workdir=$(echo $directory1 | sed -e 's/$/\/QC_frby_polish/')
#mkdir -p $workdir && cd $workdir
cd ${directory1}

CPU=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l | awk '{print $0-1}')
MEM=$(free -g | awk 'NR==2 {print $4}' | sed -e 's/\.[0-9]*//')

source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh;
module load samtools bcftools gcc;

    printf "\tRunning QC for ${reference_vcf%.vcf.gz/.numvar} \n";
    printf "\tCounting the number of changes for QV\n";
    bcftools view --threads ${CPU} -H -i 'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ov ${reference_vcf} | \
        awk -F "\t" '{print $4"\t"$5}' | \
        awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' \
            > ${reference_vcf/%.vcf.gz/.numvar};
    printf "\n ${reference_vcf}  :  #bases affected: `cat ${reference_vcf/%.vcf.gz/.numvar}`\n"


