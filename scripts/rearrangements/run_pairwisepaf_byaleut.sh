#!/usr/bin/bash

set -e

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/Pairwise/*/Align.Sebastes_aleutianus.Sebast*.PAF.all | 
    while read paf1; do
        qryspecies=$(basename "$paf1" | sed -e 's/.*Align.Sebastes_aleutianus.//' -e 's/.PAF.all$//') ;
        awk -F'\t' -v species="$qryspecies" 'BEGIN{OFS="\t"} $7>=5e+6 {$1=species"."$1; $6="Sebastes_aleutianus."$6; print $0}' $paf1 | cut -f1-12
done  >/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/Pairwise/otherrefgenomes_byaleut.chr.paf      

awk -F '\t' '$12>=40 && $7>=5e+6 {a[$6"\t"$1"\t"$5] += $10} END{for (i in a) print i"\t"a[i]}'  /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/Pairwise/otherrefgenomes_byaleut.chr.paf | sort -k2,2V -k4,4nr | awk -F'\t' '!x[$2]++' >/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/Pairwise/otherrefgenomes_byaleut.chr.paf.stranded.top

