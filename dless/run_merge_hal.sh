#!/usr/bin/sh

set -e 

[ -e TMP3 ] && rm TMP3


for type in D I V P GD GI
# NOTE	:	No Duplications; GDB; DB detected
do
	cat REF_lacvir1.HAL.bed | fgrep -w -e "$type" >1.bed
	cat REF_lacbil1.HAL.bed | fgrep -w -e "$type" >2.bed
	bedtools intersect -f 0.7 -F 0.7 -e -a 1.bed -b 2.bed -wo | cut -f1-12 | awk -F'\t' 'BEGIN{OFS="\t"} {if($3-$2 <= $9-$8) {print $1,$2,$3,$4,$5,$6} else print $7,$8,$9,$10,$11,$12}' >TMP1
	bedtools intersect -f 0.7 -F 0.7 -e -b 1.bed -a 2.bed -wo | cut -f1-12 | awk -F'\t' 'BEGIN{OFS="\t"} {if($3-$2 < $9-$8) {print $1,$2,$3,$4,$5,$6} else print $7,$8,$9,$10,$11,$12}' >>TMP1
	sort -u TMP1 | sort -t$'\t' -k1,1V -k2,2n -k3,3nr >TMP2
	rm TMP1 1.bed 2.bed
	bedtools merge -i TMP2 -c 4,5,6 -o distinct,absmin,distinct >>TMP3
	echo "Merged A and B"
	rm TMP2
	echo "$type     DONE"
done

sort -u TMP3 | sort -k1,1V -k2,2n -k3,3nr >TMP4
#rm TMP3
bedtools cluster -i TMP4 | awk -F'\t' 'BEGIN{OFS="\t"} {print $0,$3-$2}' | sort -t$'\t' -k7,7n -k8,8nr | awk -F'\t' '!x[$7]++' | cut -f1-6 >finish_HAL_SVs.bed
rm TMP4 TMP3



