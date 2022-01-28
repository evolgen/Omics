#!/usr/bin/bash

#set -e -o pipefail

outname=$4
outname_bf=${outname}.cleanstop.bf

if [ $# -ne 4 ]; then
    echo "  Need the genetic_code alignment_file.fasta out_file.fasta"
    printf "\tbash ~/RGP/scripts/selection/hyphy_cleanstop.run.sh Universal infile.fasta outfile.fasta outname \n\n"
    exit 1
fi

bash ~/RGP/scripts/selection/hyphy_cleanstop.inputcreate.sh $1 $2 $3 >${outname_bf}

hyphy ${outname_bf} >${outname}.logclean

