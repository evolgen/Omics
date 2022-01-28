#!/usr/bin/bash

set -e -o pipefail

mkdir -p submit;

#conda activate orthofinder
#module load orthofinder gcc
###/global/scratch2/rohitkolora/miniconda3/envs/orthofinder/bin/orthofinder -f $PWD/input_data -t 32 -S diamond -op 1>blast_commands.log 2>blast_commands.list
###grep '^diamond blastp -d' blast_commands.log | sed -e 's/ -p 1 / -p 20 /' >blast_commands.list

count=1
cat blast_commands.list |
#    head -n 2 |
    while read line1; do
        cat ~/RGP/scripts/orthology/dispatch_empty.sh |
            sed -e "s/name=empty$/name=${count}/" -e "s/ dispatch.log$/ D_${count}.log/" \
            >submit/${count}.sh;
        echo ${line1} >>submit/${count}.sh;
        cd submit/;
        sbatch ${count}.sh;
        cd ../;
        count=$((count+1)); sleep 5;
    done   

