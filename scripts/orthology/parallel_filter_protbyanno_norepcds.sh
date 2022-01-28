#!/usr/bin/bash -e

set -e  -o pipefail

# \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/ |
#   parallel -j 32 sh ~/RGP/scripts/orthology/parallel_filter_protbyanno.sh 

workdir=$1
#size=$2

module load samtools bedtools augustus

mkdir -p ${workdir}/FILTER;
augustus_gff="${workdir}/braker/augustus.ab_initio.gff";
species=$(echo ${workdir} | sed  -e 's/.*\/output\///' -e 's/\/.*//'); 
augustus_gff="${workdir}/braker/augustus.ab_initio.gff";
augustus_gtf="${workdir}/braker/augustus.ab_initio.gtf";
spaln_gff="${workdir}/crossprots_${species}.Genes.gff";
printf $species"\n"; 

if [ ! -e "${workdir}/FILTER/final_prots.CDS.gff" ] ; then

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

    samtools faidx ${workdir}/braker/augustus.ab_initio.aa;    
    perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${workdir}/braker/augustus.ab_initio.aa |
        grep -v -e '^[A-Z]*[A-Z]\*[A-Z]' |
        grep -B 1 -e '^M.*\*$' | 
        sed -e 's/\*$//' -e 's/\*/X/g' | 
        grep '^>' | 
        sed -e 's/^>//' | 
        sort -u \
            >${workdir}/FILTER/metstop4aug.list;
#    fgrep -w -f - ${workdir}/braker/augustus.ab_initio.aa.fai \
#    awk -F'\t' '$2>=60 {print $1}' 

    cat ${workdir}/FILTER/crossprot2aug.list ${workdir}/FILTER/metstop4aug.list |
        sort -u >${workdir}/FILTER/final_prots.list; 

    cat ${workdir}/FILTER/final_prots.list | 
        xargs samtools faidx ${workdir}/braker/augustus.ab_initio.aa |
        sed -e "/^>/ s/^>/>${species}./" >${workdir}/FILTER/final_prots.faa; 

    fgrep -w -f ${workdir}/FILTER/final_prots.list ${augustus_gff} \
        >${workdir}/FILTER/final_prots.ALL.gff; 

    fgrep -w CDS ${workdir}/FILTER/final_prots.ALL.gff | 
        awk -F'\t' -v species="${species}" 'BEGIN{OFS="\t"} {gsub("transcript_id \"","ID="species".",$NF); gsub("\".*$",";",$NF); print $0}' \
            >${workdir}/FILTER/final_prots.CDS.gff;

fi

awk '{if($6<$7){print $5"\t"$6"\t"$7"\t"$10"\t"$11"\t"9} else print $5"\t"$7"\t"$6"\t"$10"\t"$11"\t"9}' ${workdir}/REPEAT/final.genome.scf.FAS.out |
    sed -e 's/\tC$/\t-/' |
    tail -n +4 |
    sort -k1,1 -k2,2n -k3,3n \
        >${workdir}/FILTER/all.repeats.bed
bedtools merge -i ${workdir}/FILTER/all.repeats.bed >${workdir}/FILTER/all.repeats.merged.bed    

#    awk -F'\t' -v size="${size}" '$3-$2>=size' | 
bedtools intersect -v -f 0.6 -a ${workdir}/FILTER/final_prots.CDS.gff -b ${workdir}/FILTER/all.repeats.merged.bed |
    bedtools intersect -f 1 -r -u -a ${augustus_gtf} -b - |
    grep -w CDS \
        >${workdir}/FILTER/final_prots.norepeats.CDS.gtf

fgrep -v -e 'CDS' -e '^#' ${augustus_gtf} |
    bedtools intersect -u -a - -b ${workdir}/FILTER/final_prots.norepeats.CDS.gtf |
    cat - ${workdir}/FILTER/final_prots.norepeats.CDS.gtf |
    sort -k1,1 -k4,4n -k5,5nr \
        >${workdir}/FILTER/final_prots.norepeats.ALL.gtf
            
python ~/local_modules_sw/augustus/2.5.5/scripts/getAnnoFastaFromJoingenes.py -g ${workdir}/final.genome.scf.FAS \
    -o ${workdir}/FILTER/Annot_final -s True -f ${workdir}/FILTER/final_prots.norepeats.ALL.gtf

#        grep -w -e 'transcript' |
#        awk -F'\t' '{print $NF}' |
#        sed -e 's/.*=//' -e 's/;$//' -e "s/^/${species}./" |
#        sort -u \
#            >${workdir}/FILTER/final_prots.norepeats.list;

#cat ${workdir}/FILTER/final_prots.norepeats.list |
#    xargs samtools faidx ${workdir}/FILTER/final_prots.faa \
#        >${workdir}/FILTER/final_prots.norepeats.faa;

#fgrep -w -f ${workdir}/FILTER/final_prots.norepeats.list ${workdir}/FILTER/final_prots.CDS.gff \
#    >${workdir}/FILTER/final_prots.norepeats.CDS.gff;

mkdir -p ./input_filter_prots; 
cp ${workdir}/FILTER/Annot_final.aa ./input_filter_prots/${species}.faa; 
cp ${workdir}/FILTER/final_prots.norepeats.CDS.gtf ./input_filter_prots/${species}.gff;  

printf "\t$species    Done\t-\t$(grep -c '^>' ${workdir}/FILTER/Annot_final.aa)\n";



