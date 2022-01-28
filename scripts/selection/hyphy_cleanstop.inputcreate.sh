#!/usr/bin/bash

HYPHY_LIB_DIRECTORY="/global/scratch2/rohitkolora/miniconda3/envs/hyphy/lib/hyphy/"

genetic_code=$1
alignment_file=$2
out_file=$3

cat << EOF
inputRedirect = {};
inputRedirect["01"]="${genetic_code}"; // genetic code    -  Universal Vertebrate-mtDNA
inputRedirect["02"]="${alignment_file}"; // codon data      -  nexus alignment file
inputRedirect["03"]="Disallow stops"; // remove sequences with stop codons and duplicated
inputRedirect["04"]="${out_file}"; // tree            -   newick tree file

ExecuteAFile ("${HYPHY_LIB_DIRECTORY}"+"TemplateBatchFiles/CleanStopCodons.bf", inputRedirect);
EOF

#/global/scratch2/rohitkolora/miniconda3/envs/hyphy/lib/hyphy/TemplateBatchFiles/CleanStopCodons.bf

