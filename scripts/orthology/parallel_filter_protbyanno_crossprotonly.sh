#!/usr/bin/bash -e

set -e  -o pipefail

# \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/ |
#   parallel -j 32 sh ~/RGP/scripts/orthology/parallel_filter_protbyanno.sh 

workdir=$1
size=10

module load samtools bedtools augustus

mkdir -p ${workdir}/FILTER;
species=$(echo ${workdir} | sed  -e 's/.*\/output\///' -e 's/\/.*//'); 
augustus_gff="${workdir}/braker/augustus.ab_initio.gff";
augustus_gtf="${workdir}/braker/augustus.ab_initio.gtf";
spaln_gff="${workdir}/crossprots_${species}.Genes.gff";
printf $species"\n"; 

    grep -w gene ${augustus_gff} \
        >${workdir}/FILTER/augustus.ab_initio.GENE.gff;
    grep -w transcript ${augustus_gff} \
        >${workdir}/FILTER/augustus.ab_initio.TRANSCRIPT.gff;
    awk -F'\t' 'BEGIN{OFS="\t"} {if($5>$4) {print $1,$4,$5} else print $1,$5,$4}' ${spaln_gff} \
        >${workdir}/FILTER/crossprots.bed; 

    bedtools intersect -r -f 0.65 -u -a ${workdir}/FILTER/augustus.ab_initio.TRANSCRIPT.gff \
        -b ${workdir}/FILTER/crossprots.bed | 
        awk -F'\t' '{print $NF}' |
        sort -u \
            >${workdir}/FILTER/crossprot2aug.list; 

    cat ${workdir}/FILTER/crossprot2aug.list |
        xargs samtools faidx ${workdir}/braker/augustus.ab_initio.aa |
            sed -e "/^>/ s/^>/>${species}./" >${workdir}/FILTER/crossprot2aug.simpler.faa; 

    cat ${workdir}/FILTER/crossprot2aug.list | 
        fgrep -w -f - ${workdir}/braker/augustus.ab_initio.gff |
        fgrep -w -e CDS \
            >${workdir}/FILTER/crossprot2aug.simpler.CDS.gff

    awk -v species="$species" 'BEGIN{OFS="\t"} {gsub("transcript_id \"","ID="species".",$NF); gsub("\".*$",";",$NF); print $0}' \
        ${workdir}/FILTER/crossprot2aug.simpler.CDS.gff \
            >${workdir}/FILTER/crossprot2aug.simpler.CDS.edit.gff

    mkdir -p ./input_simpler_crossprots; 
    cp ${workdir}/FILTER/crossprot2aug.simpler.faa ./input_simpler_crossprots/${species}.faa; 
    cp ${workdir}/FILTER/crossprot2aug.simpler.CDS.edit.gff ./input_simpler_crossprots/${species}.gff;  

printf "\t$species    Done\t-\t$(grep -c '^>' ${workdir}/FILTER/crossprot2aug.simpler.faa)\n";



