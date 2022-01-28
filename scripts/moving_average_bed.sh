#!/usr/bin/bash

set -e

size=$1
#overlap=$2
file=$2

# awk -F'\t' '{sum+=$3} (NR%100000)==0 {print $1"\t"$2"\t"sum/100000; sum=0;}' file.bed

if [ "$#" -ne 2 ] ; then
    echo "Need three arguments : script.sh window_size #overlap_size infile.bed"
    exit 1;
fi

awk -F'\t' -v size="$size" ' {sum+=$3}
    (NR%size)==0 { print $1"\t"$2"\t"sum/size; sum=0; } ' $file

#awk -F'\t' -v size="$size" -v overlap="$overlap" '{ window = size; } (NR%window)==0 
#    { x = $3; i = NR % window; Move_Avg += (x - Z[i]) / window; Z[i] = x; print $1"\t"NR/window"\t"Move_Avg; }' $file

