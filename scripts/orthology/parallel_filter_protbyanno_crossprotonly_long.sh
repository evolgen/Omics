#!/usr/bin/bash -e

set -e  -o pipefail

# \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/ |
#   parallel -j 32 sh ~/RGP/scripts/orthology/parallel_filter_protbyanno.sh 

workdir=$1
size=10

module load samtools bedtools augustus

mkdir -p ${workdir}/FILTER;
augustus_gff="${workdir}/braker/augustus.ab_initio.gff";
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

#    cat ${workdir}/FILTER/crossprot2aug.list | 
#        sed -e "/^>/ s/^>/>${species}./" >${workdir}/FILTER/crossprot2aug.faa; 

    mkdir -p ${workdir}/FILTER/CROSSPROT/;
    cat ${workdir}/FILTER/crossprot2aug.list |
        fgrep -e -w -f - ${augustus_gff} |
        fgrep -w CDS >${workdir}/FILTER/CROSSPROT/crossprot2aug.CDS.gff

    awk '{if($6<$7){print $5"\t"$6-1"\t"$7"\t"$10"\t"$11"\t"$9} else print $5"\t"$7"\t"$6+1"\t"$10"\t"$11"\t"$9}' ${workdir}/REPEAT/final.genome.scf.FAS.out |
        sed -e 's/\tC$/\t-/' |
        tail -n +4 |
        sort -k1,1 -k2,2n -k3,3n \
            >${workdir}/FILTER/all.repeats.bed
    bedtools merge -i ${workdir}/FILTER/all.repeats.bed >${workdir}/FILTER/all.repeats.merged.bed    

    awk -F'\t' -v size="${size}" '$3-$2>=size' ${workdir}/FILTER/all.repeats.merged.bed | 
    bedtools intersect -v -f 0.6 -a ${workdir}/FILTER/CROSSPROT/crossprot2aug.CDS.gff -b - |
        bedtools intersect -f 1 -r -u -a ${augustus_gtf} -b - |
        grep -w CDS \
            >${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.CDS.gtf

    fgrep -v -e 'CDS' -e '^#' ${augustus_gtf} |
        bedtools intersect -u -a - -b ${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.CDS.gtf |
        cat - ${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.CDS.gtf |
        sort -k1,1 -k4,4n -k5,5nr \
            >${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.ALL.gtf
            
    python ~/local_modules_sw/augustus/2.5.5/scripts/getAnnoFastaFromJoingenes.py -g ${workdir}/final.genome.scf.FAS \
        -o ${workdir}/FILTER/CROSSPROT/Annot_crossprot -s True -f ${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.ALL.gtf

    sed -e "/^/ s/>/>${species}./" ${workdir}/FILTER/CROSSPROT/Annot_crossprot.aa >${workdir}/FILTER/CROSSPROT/Annot_crossprot.faa
    sed -e "/^/ s/>/>${species}./" ${workdir}/FILTER/CROSSPROT/Annot_crossprot.codingseq >${workdir}/FILTER/CROSSPROT/Annot_crossprot.fna
    awk -v species="$species" 'BEGIN{OFS="\t"} {gsub(".*transcript_id.*"g","ID="species".g",$NF); gsub("\".*$",";",$NF); print $0}' \
        ${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.CDS.gtf \
            >${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.CDS.gff

    mkdir -p ./input_filter_crossprots; 
    cp ${workdir}/FILTER/CROSSPROT/Annot_crossprot.faa ./input_filter_crossprots/${species}.faa; 
    cp ${workdir}/FILTER/CROSSPROT/crossprot2aug.norepeats.CDS.gff ./input_filter_crossprots/${species}.gff;  

printf "\t$species    Done\t-\t$(grep -c '^>' ${workdir}/FILTER/CROSSPROT/Annot_crossprot.faa)\n";



