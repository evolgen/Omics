#!/usr/bin/bash

set -e

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/codeml/sequence_files/OrthoGroup*/long/null.codeml.long.out | 
    while read nullmodel; do 
        groupname=$(dirname "$nullmodel" | sed -e 's/.*\/OrthoGroup/OrthoGroup/' -e 's/\/.*//') ;
        null_lnl=$(fgrep 'lnL' $nullmodel | sed -e 's/.*://' | awk '{print $1}') ; 
        altmodel=$(echo $nullmodel | sed -e 's/null.codeml.long.out/alt.codeml.long.out/'); 
        alt_lnl=$(fgrep 'lnL' $altmodel | sed -e 's/.*://'  | awk '{print $1}' ) ; 
        lnlval=$(echo "" | awk -v null_lnl="$null_lnl" -v alt_lnl="$alt_lnl" '{print 2*(alt_lnl - null_lnl)}') ;
        perl /global/home/users/rohitkolora/RGP/scripts/selection/calculate_chisquare_sign.pl ${lnlval} ${groupname} ; 
    
    done
