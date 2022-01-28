#!/usr/bin/bash -e

set -e  -o pipefail

# conda activate braker
# \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/ |
#   parallel -j 32 sh ~/RGP/scripts/parallel_freeze_augustus_filter.sh

workdir=$1
size=20

module load samtools bedtools augustus

if [ -e "${workdir}/referencegenome.fasta" ]; then

    mkdir -p ${workdir}/FILTERED;
    augustus_gff="${workdir}/BRAKER/augustus.ab_initio.gff";
    assembler=$(echo ${workdir} | sed  -e 's/.*\/FREEZE\///' -e 's/\/.*//')
    sample=$(echo ${workdir} | sed  -e 's/.*\/FREEZE\///' -e 's/\/$//' -e 's/.*\///')
    species="${assembler}_${sample}" 
    printf $species"\n";

    augustus_gff="${workdir}/BRAKER/augustus.ab_initio.gff";
    augustus_gtf="${workdir}/BRAKER/augustus.ab_initio.gtf";
    spaln_gff="${workdir}/crossprots_${species}.Genes.gff";

    if [ ! -e "${spaln_gff}" ]; then
    awk -F'\t' '$3=="gene"' ${workdir}/crossprots_${species}.gff >crossprots_{species}.Genes.gff
    fi

    if [ ! -e "${augustus_gtf}" ]; then
    cd ${workdir}/BRAKER/; 
    cat ${augustus_gff} | perl -ne 'if(m/\tAUGUSTUS\t/) {print $_;}' | perl /global/home/users/rohitkolora/local_modules_sw/augustus/2.5.5/scripts/gtf2gff.pl --printExon --out=${augustus_gtf};
    fi
    
    if [ ! -e "${workdir}/FILTERED/final_prots.CDS.gff" ] ; then

    grep -w gene ${augustus_gff} \
        >${workdir}/FILTERED/augustus.ab_initio.GENE.gff;
    grep -w transcript ${augustus_gff} \
        >${workdir}/FILTERED/augustus.ab_initio.TRANSCRIPT.gff;
    awk -F'\t' 'BEGIN{OFS="\t"} {if($5>$4) {print $1,$4,$5} else print $1,$5,$4}' ${spaln_gff} \
        >${workdir}/FILTERED/crossprots.bed; 

    bedtools intersect -r -f 0.65 -u -a ${workdir}/FILTERED/augustus.ab_initio.TRANSCRIPT.gff \
        -b ${workdir}/FILTERED/crossprots.bed | 
        awk -F'\t' '{print $NF}' |
        sort -u \
            >${workdir}/FILTERED/crossprot2aug.list; 

    if [ ! -e "${workdir}/BRAKER/augustus.ab_initio.aa" ]; then
        samtools faidx ${workdir}/BRAKER/augustus.ab_initio.aa;    
    fi    
    perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${workdir}/BRAKER/augustus.ab_initio.aa |
        grep -v -e '^[A-Z]*[A-Z]\*[A-Z]' |
        grep -B 1 -e '^M.*\*$' | 
        sed -e 's/\*$//' -e 's/\*/X/g' | 
        grep '^>' | 
        sed -e 's/^>//' | 
        sort -u \
            >${workdir}/FILTERED/metstop4aug.list;

    cat ${workdir}/FILTERED/crossprot2aug.list ${workdir}/FILTERED/metstop4aug.list |
        sort -u >${workdir}/FILTERED/final_prots.list; 

    cat ${workdir}/FILTERED/final_prots.list | 
        xargs samtools faidx ${workdir}/BRAKER/augustus.ab_initio.aa |
        sed -e "/^>/ s/^>/>${species}./" >${workdir}/FILTERED/final_prots.faa; 
    samtools faidx ${workdir}/FILTERED/final_prots.faa;

    fgrep -w -f ${workdir}/FILTERED/final_prots.list ${augustus_gff} \
        >${workdir}/FILTERED/final_prots.ALL.gff; 

    fgrep -w CDS ${workdir}/FILTERED/final_prots.ALL.gff | 
        awk -F'\t' -v species="${species}" 'BEGIN{OFS="\t"} {gsub("transcript_id \"","ID="species".",$NF); gsub("\".*$",";",$NF); print $0}' \
            >${workdir}/FILTERED/final_prots.CDS.gff;

    fi

    awk 'BEGIN{OFS="\t"} {if($6<$7){print $5,$6-1,$7,$10,$11,$9} else print $5,$7,$6+1,$10,$11,$9}' ${workdir}/REPEAT/referencegenome.fasta.out |
    sed -e 's/\tC$/\t-/' |
    tail -n +4 |
    sort -k1,1 -k2,2n -k3,3n \
        >${workdir}/FILTERED/all.repeats.bed;

    bedtools merge -i ${workdir}/FILTERED/all.repeats.bed >${workdir}/FILTERED/all.repeats.merged.bed    

    awk -F'\t' -v size="${size}" '$3-$2>=size' ${workdir}/FILTERED/all.repeats.merged.bed | 
    bedtools intersect -v -f 0.3 -a ${workdir}/FILTERED/final_prots.CDS.gff -b - |
        bedtools intersect -f 1 -r -u -a ${augustus_gtf} -b - |
        grep -w CDS \
            >${workdir}/FILTERED/final_prots.norepeats.CDS.gtf

    fgrep -v -e 'CDS' -e '^#' ${augustus_gtf} |
    bedtools intersect -u -a - -b ${workdir}/FILTERED/final_prots.norepeats.CDS.gtf |
    cat - ${workdir}/FILTERED/final_prots.norepeats.CDS.gtf |
    sort -k1,1 -k4,4n -k5,5nr \
        >${workdir}/FILTERED/final_prots.norepeats.ALL.gtf
            
    python ~/local_modules_sw/augustus/2.5.5/scripts/getAnnoFastaFromJoingenes.py -g ${workdir}/referencegenome.fasta \
        -o ${workdir}/FILTERED/Annot_final -s True -f ${workdir}/FILTERED/final_prots.norepeats.ALL.gtf;
    samtools faidx ${workdir}/FILTERED/Annot_final.codingseq && ${workdir}/FILTERED/Annot_final.aa;

fi

if [ ! -e "${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gff" ]; then
    if [ -e "${workdir}/referencegenome.fasta" ]; then
    mkdir -p ${workdir}/FILTERED/CROSSPROT/ && cd ${workdir}/FILTERED/CROSSPROT/;
    cat ${workdir}/FILTERED/crossprot2aug.list |
        fgrep -e -w -f - ${augustus_gff} |
        fgrep -w CDS >${workdir}/FILTERED/CROSSPROT/crossprot2aug.CDS.gff

    awk '{if($6<$7){print $5"\t"$6-1"\t"$7"\t"$10"\t"$11"\t"$9} else print $5"\t"$7"\t"$6+1"\t"$10"\t"$11"\t"$9}' ${workdir}/REPEAT/referencegenome.fasta.out |
        sed -e 's/\tC$/\t-/' |
        tail -n +4 |
        sort -k1,1 -k2,2n -k3,3n \
            >${workdir}/FILTERED/all.repeats.bed
    bedtools merge -i ${workdir}/FILTERED/all.repeats.bed >${workdir}/FILTERED/all.repeats.merged.bed

    awk -F'\t' -v size="${size}" '$3-$2>=size' ${workdir}/FILTERED/all.repeats.merged.bed |
        bedtools intersect -v -f 0.6 -a ${workdir}/FILTERED/CROSSPROT/crossprot2aug.CDS.gff -b - |
        bedtools intersect -f 1 -r -u -a ${augustus_gtf} -b - |
        grep -w CDS \
            >${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gtf

    fgrep -v -e 'CDS' -e '^#' ${augustus_gtf} |
        bedtools intersect -u -a - -b ${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gtf |
        cat - ${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gtf |
        sort -k1,1 -k4,4n -k5,5nr \
            >${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.ALL.gtf

    python ~/local_modules_sw/augustus/2.5.5/scripts/getAnnoFastaFromJoingenes.py -g ${workdir}/referencegenome.fasta \
        -o ${workdir}/FILTERED/CROSSPROT/Annot_crossprot -s True -f ${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.ALL.gtf;
    samtools faidx ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.codingseq && ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.aa;

    sed -e "/^/ s/>/>${species}./" ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.aa >${workdir}/FILTERED/CROSSPROT/Annot_crossprot.faa
    sed -e "/^/ s/>/>${species}./" ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.codingseq >${workdir}/FILTERED/CROSSPROT/Annot_crossprot.fna
    samtools faidx ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.faa && samtools faidx ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.codingseq
    awk -v species="$species" 'BEGIN{OFS="\t"} {gsub(".*transcript_id.*\"g","ID="species".g",$NF); gsub("\".*$",";",$NF); print $0}' \
        ${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gtf \
            >${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gff;

    awk -F'\t' species="${species}" 'gsub("\";.*",";",$0); gsub("transcript_id.*\"g","ID="species".g",$0); print $0' ${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gff >${workdir}/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.GFF;
    sed -e "/^>/ s/>/>${species}./" ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.faa >${workdir}/FILTERED/CROSSPROT/Annot_crossprot.FAA;

printf "\t$species    Done\t-\t$(grep -c '^>' ${workdir}/FILTERED/CROSSPROT/Annot_crossprot.faa)\n";
    fi
fi

