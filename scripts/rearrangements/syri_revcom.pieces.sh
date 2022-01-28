#!/usr/bin/bash

set -e -o pipefail

module load samtools

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/WTDBG/Sebastolobus_alascanus/long_seq.fa.orderby.Sebastes_aleutianus"
printf "  Processing $file\n";
samtools faidx $file;
printf "" >${file}.listedbyaleu;
for num in 3; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf "\tReverse complementing  :  ${chromo}\t${num}\n";
    printf ">${chromo}\n" >>${file}.listedbyaleu;
    samtools faidx ${file} $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.listedbyaleu;
done
samtools faidx ${file}.listedbyaleu;
cut -f1 ${file}.listedbyaleu.fai >${file}.listedbyaleu.names;
cut -f1 ${file/%.orderby.Sebastes_aleutianus/}.fai | 
    grep -v -f ${file}.listedbyaleu.names | 
    xargs samtools faidx ${file/%.orderby.Sebastes_aleutianus/} >>${file}.listedbyaleu; 
samtools faidx ${file}.listedbyaleu;        
printf "\tCreating the new edited final file\n";        
cut -f1 ${file}.fai | 
    xargs samtools faidx ${file}.listedbyaleu >${file}.edit;
samtools faidx ${file}.edit;
rm -f ${file}.listedbyaleu.names ${file}.listedbyaleu ${file}.listedbyaleu.fai;
echo "Done with    ${file}.edit";


file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/WTDBG/Sebastes_entomelas/long_seq.fa.orderby.Sebastes_aleutianus"
printf "  Processing $file\n";
samtools faidx $file;
printf "" >${file}.listedbyaleu;
for num in 1 10 11 3 4 6 9 ; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf "\tReverse complementing  :  ${chromo}\t${num}\n";
    printf ">${chromo}\n" >>${file}.listedbyaleu;
    samtools faidx ${file} $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.listedbyaleu;
done
samtools faidx ${file}.listedbyaleu;
cut -f1 ${file}.listedbyaleu.fai >${file}.listedbyaleu.names;
cut -f1 ${file/%.orderby.Sebastes_aleutianus/}.fai | 
    grep -v -f ${file}.listedbyaleu.names | 
    xargs samtools faidx ${file/%.orderby.Sebastes_aleutianus/} >>${file}.listedbyaleu; 
samtools faidx ${file}.listedbyaleu;        
printf "\tCreating the new edited final file\n";        
cut -f1 ${file}.fai | 
    xargs samtools faidx ${file}.listedbyaleu >${file}.edit;
samtools faidx ${file}.edit;
rm -f ${file}.listedbyaleu.names ${file}.listedbyaleu ${file}.listedbyaleu.fai;
echo "Done with    ${file}.edit";

file="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_miniatus/long_seq.fa.orderby.Sebastes_aleutianus"
printf "  Processing $file\n";
samtools faidx $file;
printf "" >${file}.listedbyaleu;
for num in 10 6 8 9 12 14 2 ; do
    chromo=$(awk -v num="$num" 'NR==num' ${file}.fai | cut -f1);
    printf "\tReverse complementing  :  ${chromo}\t${num}\n";
    printf ">${chromo}\n" >>${file}.listedbyaleu;
    samtools faidx ${file} $chromo |
        grep -v '^>' | tr "\n" " " | sed -e 's/ //g' |
        tr "[ATGCatgc]" "[TACGtacg]" | rev >>${file}.listedbyaleu;
done
samtools faidx ${file}.listedbyaleu;
cut -f1 ${file}.listedbyaleu.fai >${file}.listedbyaleu.names;
cut -f1 ${file/%.orderby.Sebastes_aleutianus/}.fai | 
    grep -v -f ${file}.listedbyaleu.names | 
    xargs samtools faidx ${file/%.orderby.Sebastes_aleutianus/} >>${file}.listedbyaleu; 
samtools faidx ${file}.listedbyaleu;        
printf "\tCreating the new edited final file\n";        
cut -f1 ${file}.fai | 
    xargs samtools faidx ${file}.listedbyaleu >${file}.edit; 
samtools faidx ${file}.edit;
rm -f ${file}.listedbyaleu.names ${file}.listedbyaleu ${file}.listedbyaleu.fai;
echo "Done with    ${file}.edit";


