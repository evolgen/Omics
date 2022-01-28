#!/usr/bin/bash

set -e

present_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge";

count=0; 
cat /global/scratch2/rohitkolora/Rockfish/Data/Greg/busco/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.ufboot | 
    sort | uniq -c | 
    awk '{print $2"\t"$1}' | while read line; do 
        count=$((count+1)); 
        mkdir -p ${present_dir}/tree_multi/${count}; 
        tree=$(echo $line | sed -e 's/ .*//'); 
        weight=$(echo $line | sed -e 's/.* //'); 
        echo $tree >${present_dir}/tree_multi/${count}/tree.nwk; 
        echo $weight >${present_dir}/tree_multi/${count}/weight.txt; 
        echo $count ;
done


