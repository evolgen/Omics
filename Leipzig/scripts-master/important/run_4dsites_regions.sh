#!/usr/bin/sh

set -e

#   sh run_4dsites_regions.sh /scr/bloodymary/rohit/Lacerta_viridis/Large_SVs/Lacerta_viridis_collinear.bed /scr/bloodymary/rohit/Lacerta_viridis/Selection/Divergence_time_MA.txt 

file1=$1
file2=$2

fileout=$(echo $file1 | sed -e 's/.*\///' -e 's/.bed/_4d-divergence.txt/')
fileout_average=$(echo $file1 | sed -e 's/.*\///' -e 's/.bed/_4d-divergence.average.txt/')
fileout_median=$(echo $file1 | sed -e 's/.*\///' -e 's/.bed/_4d-divergence.median.txt/')


cat $file1 | while read -r line; do echo $line | sed -e 's/ /\t/g' | tr "\n" "\t"; echo $line | sed -e 's/ /\t/g' | bedtools intersect -u -a /scr/bloodymary/rohit/Lacerta_viridis/annotation/Lacerta_viridis_cdsbasedannotation.bed_names -b - | cut -f5 | while read -r transcript; do printf $(echo $transcript | fgrep -w -f - $file2 | cut -f7)"\n"; done; echo; done | sed -e '/^$/d' | awk 'NF!=6' | tr "\n" "\t" | sed -e 's/Lvir_/\nLvir_/g' | sed '/^$/d' | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,NF-6}' >tmp1

cat $file1 | while read -r line; do echo $line | sed -e 's/ /\t/g' | tr "\n" "\t"; echo $line | sed -e 's/ /\t/g' | bedtools intersect -u -a /scr/bloodymary/rohit/Lacerta_viridis/annotation/Lacerta_viridis_cdsbasedannotation.bed_names -b - | cut -f5 | while read -r transcript; do printf $(echo $transcript | fgrep -w -f - $file2 | cut -f7)"\n"; done; echo; done | sed -e '/^$/d' | awk 'NF!=6' | tr "\n" "\t" | sed -e 's/Lvir_/\nLvir_/g' | sed '/^$/d' | cut -f7- | while read line2; do echo $line2 | sed -e 's/ /\n/g' | awk -F'\t' '{sum+=$0} END {print sum}'; done >tmp2

cat $file1 | while read -r line; do echo $line | sed -e 's/ /\t/g' | tr "\n" "\t"; echo $line | sed -e 's/ /\t/g' | bedtools intersect -u -a /scr/bloodymary/rohit/Lacerta_viridis/annotation/Lacerta_viridis_cdsbasedannotation.bed_names -b - | cut -f5 | while read -r transcript; do printf $(echo $transcript | fgrep -w -f - $file2 | cut -f7)"\n"; done; echo; done | sed -e '/^$/d' | awk 'NF!=6' | tr "\n" "\t" | sed -e 's/Lvir_/\nLvir_/g' | sed '/^$/d' | cut -f7- | while read line2; do echo $line2 | sed -e 's/ /\n/g' | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }' ; done >tmp3

paste tmp1 tmp2 | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$8/$7}' >$fileout_average
paste tmp1 tmp3 | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$8/$7}' >$fileout_median


rm tmp1 tmp2 tmp3


