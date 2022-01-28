#!/usr/bin/bash

set -e -o pipefail

file_tsv=$@

name=$(dirname $file_tsv | sed -e 's/.*searches\///' -e 's/\//__/')
file_peaks=$(echo $file_tsv | sed -e 's/.tsv$/.51.cvg.peaks/');
file_count=$(echo $file_tsv | sed -e 's/.tsv$/.count.txt/');

n_chr=24
pl_length=6
read_length=150

if [[ -e "$file_peaks" && -e "$file_count" ]]; then

base_coverage=$(cat $file_peaks | grep '^#haploid_fold_coverage' | awk '{print $2}');
tel_coverage=$(awk '{print $1}' $file_count);

rel_coverage=$((tel_coverage/base_coverage));
#printf "rel_coverage : ${rel_coverage}";

telomore_length=$(awk -v rel_coverage="$rel_coverage" -v n_chr="$n_chr" -v pl_length="$pl_length" -v read_length="$read_length" -v base_coverage="$base_coverage" '{print (rel_coverage*((read_length+pl_length-1)/(2*n_chr))) / 1e+6 }' ${file_count});

printf "$name\t${telomore_length}\n";

fi

