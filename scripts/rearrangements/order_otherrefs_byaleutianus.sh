#!/usr/bin/bash

set -e -o pipefail

if [ "$#" -ne 2 ]; then
    echo "   Need an input paf file and output file";
    exit 1;
fi

columncheck=$(awk -F'\t' '{print NF}' $file1 | sort -u);
columncount=$(awk -F'\t' '{print NF}' $file1 | sort -u | wc -l);

if [ "$columncheck" -ne 12 ] || [ "$columncount" -ne 1 ] ; then
    echo "   Strictly need a 12-column paf file validity";
    echo "   You have $columncount column types and column numbers = $columncheck";
    exit 1;
fi

work_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE"

\ls -d ${work_dir}/FALCON/sm2}/Pairwise/{asm2}/align.{asm1}.{sp2}.paf2

