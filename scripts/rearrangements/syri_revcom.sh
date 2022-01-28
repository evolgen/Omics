#!/usr/bin/bash

set -e -o pipefail

module load samtools
file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/long_seq.fa"
samtools faidx $file;
cp $file ${file}.edit;
echo "Done with    $file"

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_miniatus/long_seq.fa.orderby.Sebastes_aleutianus"
samtools faidx $file;
printf "" >${file}.edit;
for num in 10 11 17 18 19 2 22 4 7; do 
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1); 
    printf ">${chromo}\n" >>${file}.edit; 
    samtools faidx $file $chromo | 
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' | 
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.edit; 
done
samtools faidx ${file}.edit;
cut -f 1 ${file}.edit.fai | grep -v -f - ${file}.fai | cut -f1 | xargs samtools faidx ${file} >>${file}.edit
samtools faidx ${file}.edit;
echo "Done with    $file"

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_pinniger/long_seq.fa.orderby.Sebastes_aleutianus"
samtools faidx $file;
printf "" >${file}.edit;
for num in 11 12 14 22 24 3 6 8; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf ">${chromo}\n" >>${file}.edit;
    samtools faidx $file $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.edit;
done
samtools faidx ${file}.edit;
cut -f 1 ${file}.edit.fai | grep -v -f - ${file}.fai | cut -f1 | xargs samtools faidx ${file} >>${file}.edit
samtools faidx ${file}.edit;
echo "Done with    $file"

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_rosaceus/long_seq.fa.orderby.Sebastes_aleutianus"
samtools faidx $file;
printf "" >${file}.edit;
for num in 11 12 15 17 19 22 23 6 ; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf ">${chromo}\n" >>${file}.edit;
    samtools faidx $file $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.edit;
done
samtools faidx ${file}.edit;
cut -f 1 ${file}.edit.fai | grep -v -f - ${file}.fai | cut -f1 | xargs samtools faidx ${file} >>${file}.edit
samtools faidx ${file}.edit;
echo "Done with    $file"

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_umbrosus/long_seq.fa.orderby.Sebastes_aleutianus"
samtools faidx $file;
printf "" >${file}.edit;
for num in 10 14 15 17 18 4 7 8; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf ">${chromo}\n" >>${file}.edit;
    samtools faidx $file $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.edit;
done
samtools faidx ${file}.edit;
cut -f 1 ${file}.edit.fai | grep -v -f - ${file}.fai | cut -f1 | xargs samtools faidx ${file} >>${file}.edit
samtools faidx ${file}.edit;
echo "Done with    $file"

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/WTDBG/Sebastes_entomelas/long_seq.fa.orderby.Sebastes_aleutianus"
samtools faidx $file;
printf "" >${file}.edit;
for num in 10 12 13 14 15 17 18 19 20 24 5 8 ; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf ">${chromo}\n" >>${file}.edit;
    samtools faidx $file $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.edit;
done
samtools faidx ${file}.edit;
cut -f 1 ${file}.edit.fai | grep -v -f - ${file}.fai | cut -f1 | xargs samtools faidx ${file} >>${file}.edit
samtools faidx ${file}.edit;
echo "Done with    $file"

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/WTDBG/Sebastolobus_alascanus/long_seq.fa.orderby.Sebastes_aleutianus"
samtools faidx $file;
printf "" >${file}.edit;
for num in 1 3 ; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf ">${chromo}\n" >>${file}.edit;
    samtools faidx $file $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.edit;
done
samtools faidx ${file}.edit;
cut -f 1 ${file}.edit.fai | grep -v -f - ${file}.fai | cut -f1 | xargs samtools faidx ${file} >>${file}.edit
samtools faidx ${file}.edit;
echo "Done with    $file"


