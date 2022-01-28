#!/usr/bin/bash

set -e

if [ "$#" -ne 2 ]; then
    echo " Use program as : bash split_files_maf.sh infile.maf 10"
    exit 1 ;
fi

printf "\nCreating folder and removing files\n" ;
mkdir -p maf_now_files1 ;
find ./maf_now_files1/ -type f -name "*.maf" | xargs rm -f ;
find ./maf_now_files1/ -type f -name "*.lst" | xargs rm -f ;
find ./maf_now_files1/ -type f -name "*.txt" | xargs rm -f ;
find ./maf_now_files1/ -type f -name "*.info" | xargs rm -f ;

infile_maf1=$1
num_parts1=$2
num_parts1=$((num_parts1-1))

uniq_id1=$(shuf -i 1-1000000 -n 1)
total_lines1=$(wc -l $infile_maf1)

fgrep -n 'a score=' $infile_maf1 | cut -f1 -d':' >./maf_now_files1/${uniq_id1}.lst
total_blocks1=$(wc -l ./maf_now_files1/${uniq_id1}.lst | sed -e 's/ .*//')

size_factor1=$(echo "" | awk -v num_parts="$num_parts1" -v total_blocks="$total_blocks1" '{print int(total_blocks/num_parts)}')

printf "  CustomID=${uniq_id1}\tBlocks=${total_blocks1}\tSizes=${size_factor1}\n" ;

awk -v size_factor1="$size_factor1" 'NR % size_factor1 == 0 {print }' ./maf_now_files1/${uniq_id1}.lst | sed -e '/^$/d' >./maf_now_files1/${uniq_id1}_line_numbers1.txt ;
printf "  All line numbers extracted\n" ;

start1=1
idx1=1

num1=$(head -n 1 ./maf_now_files1/${uniq_id1}_line_numbers1.txt) ; # first_maf1
printf "  MAF countdown ::: ${idx1} " ;
head -n $((num1-1)) ${infile_maf1} >./maf_now_files1/$idx1.maf ;
start1=$(( num1 )) ;

printf "Index : Start : End \n" >./maf_now_files1/${uniq_id1}.info;
printf "${idx1} : $num1 : $start1 \n" >>./maf_now_files1/${uniq_id1}.info;
wc -l ./maf_now_files1/${idx1}.maf >>./maf_now_files1/${uniq_id1}.info;

tail -n +2 ./maf_now_files1/${uniq_id1}_line_numbers1.txt | while read num1;  # for num1 in "${line_numbers1[@]}";
do
    idx1=$(( idx1 + 1 )) ;
    printf "##maf version=1\n\n" >./maf_now_files1/$idx1.maf ; 
    tail -n +$(( start1-1 )) ${infile_maf1} | head -n $(( num1-start1 )) >>./maf_now_files1/${idx1}.maf ;
    printf "${idx1} : $num1 : $start1 \n" >>./maf_now_files1/${uniq_id1}.info;
    printf "${idx1} " ;
    wc -l ./maf_now_files1/${idx1}.maf >>./maf_now_files1/${uniq_id1}.info;
    start1=$(( num1 + 1 )) ;
done ;
idx1=$(( num_parts1 + 1 )) ;
printf "${idx1}\n" ;

num1=$(tail -n 1 ./maf_now_files1/${uniq_id1}_line_numbers1.txt) ; 
printf "##maf version=1\n\n" >./maf_now_files1/${idx1}.maf ;
tail -n +${num1} ${infile_maf1} >>./maf_now_files1/${idx1}.maf ;
wc -l ./maf_now_files1/${idx1}.maf >>./maf_now_files1/${uniq_id1}.info;
printf "Finished creating ${idx1} MAFs\n\n"

rm -f ./maf_now_files1/${uniq_id1}.lst ./maf_now_files1/${uniq_id1}_line_numbers1.txt ./maf_now_files1/${uniq_id1}.info ;

