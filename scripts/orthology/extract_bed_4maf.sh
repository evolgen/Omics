#!/usr/bin/bash

set -e

if [ "$#" -ne 1 ] ; then
    echo "Use as script.sh maffile " ;
    exit 1
fi    

maffile=$1
bedfile=$(echo $maffile | sed -e 's/.maf$/.bed/' -e 's/.MAF$/.bed/' -e 's/.Maf$/.bed/')
printf "$bedfile " 

awk '$1=="s"' ${maffile} | awk '{OFS="\t"}{if ($5=="-") {print $2,$3-$4,$3+1,".\t.",$5} else print $2,$3,$3+$4+1,".\t.",$5}' >${bedfile} ;

