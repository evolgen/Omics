#!/usr/bin/bash

HYPHY_LIB_DIRECTORY="/global/scratch2/rohitkolora/miniconda3/envs/hyphy/lib/hyphy/TemplateBatchFiles/"

genetic_code=$1
alignment_file=$2
tree_file=$3
outfile=$4

cat << EOF
inputRedirect = {};
inputRedirect["01"]="${alignment_file}"; // genetic code    -  Universal Vertebrate-mtDNA
inputRedirect["02"]="${genetic_code}"; // codon data      -  nexus alignment file
inputRedirect["03"]="${tree_file}"; // tree            -   newick tree file
inputRedirect["04"]="Codon Model"; // 
inputRedirect["05"]="MG94xREV"; // 
inputRedirect["06"]="Run Custom"; // 1 = All models complete selection
inputRedirect["07"]="Proportional"; // (1):[Constant] (2):[Proportional] (3):[Nonsynonymous] (4):[Dual] (5):[Lineage Dual]
inputRedirect["08"]="";
inputRedirect["09"]="Independent Discrete";
inputRedirect["10"]="Default"; // default initial values for rate distribution
inputRedirect["11"]="3"; // #synonymous rate classes
inputRedirect["12"]="${outfile}.summary";
inputRedirect["13"]="";

ExecuteAFile ("${HYPHY_LIB_DIRECTORY}"+"dNdSRateAnalysis.bf", inputRedirect);
EOF

