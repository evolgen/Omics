#!/usr/bin/bash

set -e -o pipefail

###   sh ~/RGP/scripts/orthology/process_ovlexons.sh \
###       /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/{species}/{smp}/assembly/masurca/FILTERED/CROSSPROT/FINALE/ovl.norepeats.CDS.gtf \
###       /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/{species}/{smp}/assembly/masurca/FILTERED/CROSSPROT/crossprot2aug.norepeats.CDS.gtf

module load bedtools samtools;
ovl_cdsgtf=$@
workdir=$(dirname $ovl_cdsgtf);
norepcds_gtf=${workdir}/../crossprot2aug.norepeats.CDS.gtf;
asmb_edit=$(dirname $ovl_cdsgtf | sed -e 's/\/assembly\/masurca\/FILTERED\/.*/\/assembly\/masurca\/final.genome.scf.FAS/')

cd ${workdir};
printf " creating separate files for each overlap\n"
awk -F'\t' '{ if (!($NF in a)) a[$NF] = $0; } END { for (i in a) print a[i]}' \
    ${ovl_cdsgtf} \
    >${ovl_cdsgtf}.sing
awk -F'\t' '{if (x[$NF]) { x_count[$NF]++; print $0; if (x_count[$NF] == 1) { print x[$NF] } } x[$NF] = $0}' \
    ${ovl_cdsgtf} | \
    bedtools intersect -wo -b /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/frozen.all.gene.gtf -a - | \
    awk -F'\t' '{if (x[$13$14$(NF-1)]) { x_count[$13$14$(NF-1)]++; print $0; if (x_count[$13$14$(NF-1)] == 1) { print x[$13$14$(NF-1)] } } x[$13$14$(NF-1)] = $0}' | \
    cut -f1-13 \
    >${ovl_cdsgtf}.reps;

files_ovl=$(awk -F'\t' '{print $NF}' ${ovl_cdsgtf}.reps | sort -u | wc -l)
count_ovl=$(cat ${ovl_cdsgtf}.reps | wc -l)

if [[ "${count_ovl}" == "0" ]]; then
    echo "" >NO-OVERLAPS.ERR
    printf "\tSorry could not find gene info based overlaps\n";
    exit 1;
fi    

printf " Total overlaps = ${count_ovl}\n Total files = ${files_ovl}\n ";

mkdir -p tempo && cd tempo && rm -f ./*;
printf "" >multi.cds.transcript;

printf "  Creating the bedfiles\n";
count=1;
awk -F'\t' '{ print >> $NF".bed12" ; close($NF".bed12")}' ${ovl_cdsgtf}.reps;
\ls -d ./[0-9]*.bed12 | sort -V | #head -n 30 |  
    while read bedfile1; do
        gtffile2=$(basename ${bedfile1} | sed -e 's/.bed12$/.gtf/');
        name1=$(basename ${bedfile1} | sed -e 's/.bed12$//' -e 's/^/edG/');
        printf " $count";
        cut -f4 $bedfile1 | \
            fgrep -w -f - ${norepcds_gtf} | \
            sed -e "s/\ttranscript_id .*/\ttranscript_id \"${name1}.t1\"; gene_id \"${name1}\";/" \
                >${gtffile2}.cds;
            sed -e 's/\tCDS\t/\texon\t/' ${gtffile2}.cds >${gtffile2}.exon;
            cat ${gtffile2}.cds ${gtffile2}.exon | sort -k1,1V -k4,4n -k3,3V -k5,5n >${gtffile2}.1;
            count_chr=$(cut -f1 ${gtffile2}.1 | sort -u | wc -l);
            if [[ "${count_chr}" == 1 ]]; then
                strand=$(head -n 1 ${gtffile2}.1 | cut -f7);
                chr=$(head -n 1 ${gtffile2}.1 | cut -f1);
                if [[ "$strand" == "+" ]]; then
                    start=$(head -n 1 ${gtffile2}.1 | cut -f4);
                    start_end=$((start+2));
                    stop=$(tail -n 1 ${gtffile2}.1 | cut -f5);
                    stop_begin=$((stop-2));
                    printf "${chr}\tAUGUSTUS\tgene\t${start}\t${stop}\t0.01\t${strand}\t.\t${name1}\n" >${gtffile2}.2;
                    printf "${chr}\tAUGUSTUS\ttranscript\t${start}\t${stop}\t0.01\t${strand}\t.\t${name1}\n" >>${gtffile2}.2;
                    printf "${chr}\tAUGUSTUS\tstart_codon\t${start}\t${start_end}\t0.01\t${strand}\t.\ttranscript_id \"${name1}.t1\"; gene_id \"${name1}\";\n" >>${gtffile2}.2;
                    cat ${gtffile2}.1 >>${gtffile2}.2;
                    printf "${chr}\tAUGUSTUS\tstop_codon\t${stop_begin}\t${stop}\t0.01\t${strand}\t.\ttranscript_id \"${name1}.t1\"; gene_id \"${name1}\";\n" >>${gtffile2}.2;
                else
                    start=$(head -n 1 ${gtffile2}.1 | cut -f4);
                    start_end=$((start+2));
                    stop=$(tail -n 1 ${gtffile2}.1 | cut -f5);
                    stop_begin=$((stop-2));
                    printf "${chr}\tAUGUSTUS\tgene\t${start}\t${stop}\t0.01\t${strand}\t.\t${name1}\n" >${gtffile2}.2;
                    printf "${chr}\tAUGUSTUS\ttranscript\t${start}\t${stop}\t0.01\t${strand}\t.\t${name1}\n" >>${gtffile2}.2;
                    printf "${chr}\tAUGUSTUS\tstop_codon\t${start}\t${start_end}\t0.01\t${strand}\t.\ttranscript_id \"${name1}.t1\"; gene_id \"${name1}\";\n" >>${gtffile2}.2;
                    cat ${gtffile2}.1 >>${gtffile2}.2;
                    printf "${chr}\tAUGUSTUS\tstart_codon\t${stop_start}\t${stop}\t0.01\t${strand}\t.\ttranscript_id \"${name1}.t1\"; gene_id \"${name1}\";\n" >>${gtffile2}.2;
                fi
            else
                cat ${bedfile1} | cut -f4 | \
                fgrep -w -f - ${norepcds_gtf} \
                    >${gtffile2}.3;
                awk '{ print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$7 }' ${gtffile2}.3 \
                    >${gtffile2}.3.bed;
                printf ">${name1}.t1\n" >>multi.cds.transcript;
                bedtools getfasta -s -fi ${asmb_edit} \
                    -bed ${gtffile2}.3.bed \
                    | grep -v '^>' | tr "\n" " " | sed -e 's/ //g' >>multi.cds.transcript;
                echo "" >>multi.cds.transcript;
            fi
            count=$((count+1));
    done; echo; 
    cp ${workdir}/tempo/multi.cds.transcript ${workdir}/ && cd ${workdir}/;

printf "\n\tDone taking things out - now go annotate the artificial transcripts\n";
cat ./tempo/*.gtf.2 ./tempo/*.gtf.3 >${ovl_cdsgtf}.ext;
 
