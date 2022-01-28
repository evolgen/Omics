#!/usr/bin/bash

set -e

filecpgreport=$1

echo $filecpgreport ;

fgrep -e 'ID ' -e 'CpG island' -e '/Percent CG=' -e '/ObsExp=' $filecpgreport | 
    sed -e 's/^ID   /ID::/' -e '/^ID/ s/ .*//' -e '/^FT / s/.*=//' -e '/^FT / s/.* //' |
    awk '{if($1 ~ /\.\./) {print "pos::"$0} else print $0} ' | 
    sed -e 's/\.\./\t/' | 
    awk 'ORS = $1 ~ /^[0-9p]/ ? "\t" : "\n"' | 
    sed -e 's/ID::/\nID::/g' -e 's/\tpos::/\npos::/g' -e 's/\t$/\n/' -e '/^$/d' >${filecpgreport}.extract

perl ~/RGP/scripts/cpgreport_to_bed.pl ${filecpgreport}.extract |
    awk -F'\t' 'BEGIN{OFS="\t"} {$2=$2-1; print $0}' >${filecpgreport}.bed   



