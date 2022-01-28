#!/usr/bin/bash

#set -e -o pipefail

cds=$1
extract=$2
gen_code=$3

if [ "$#" -ne "3" ]; then
    echo "  Need three arguments here : script.sh infile outfile geneticcode " ;
    exit 1;
fi

namer=$(date '+%d%m%Y_%H.%M.%S') ;

cd $(dirname "${extract}") ;
#cp $cds ${namer}.fas ;

#printf "\tRunning pre hyphy\n";
#rm -f error.log prealign.log ;
#hyphy /global/scratch2/rohitkolora/Software/hyphy-develop/hyphy-analyses-master/codon-msa/pre-msa.bf --input ${namer}.fas --code ${gen_code} --E 0.9 1>prealign.log 2>&1 ;

#mkdir -p tester && mv ${namer}.fas ${namer}.fas_nuc.fas ${namer}.fas_protein.fas tester ;
#printf "\tPre hyphy'ing\n";

rm -f ${extract}.orfout ;

#infile="$(echo $(dirname "${extract}") | sed -e 's/$/\/errors.log/')" ;
##if [[ ! -f "$infile" ]]; then
##   printf "\tNo ORF finding\n";
##    bash ~/RGP/scripts/selection/extract_filter4cds_4snakemake.noorf.sh ${cds} ${extract} ;
##elif [[ -f "$infile" ]]; then
    printf "\tORF finder used\n";
    bash ~/RGP/scripts/selection/extract_filter4cds_4snakemake.20200622.sh ${cds} ${extract} ;
##fi    

printf "\tReady for getting Hyphy'ed\n";

