#!/usr/bin/bash

set -e -o pipefail

if [ "$#" -ne 2 ]; then
    echo "    Need One input & One output fasta file"
    exit 1;
fi

infasta=$1
outfasta=$2

awk '{for(x=1;x<=NF;x++)if($x~/^>.*/){sub(/^>.*/,">Chr"++i)}}1' $infasta >$outfasta; 
printf "\tProcessed $1\n";


