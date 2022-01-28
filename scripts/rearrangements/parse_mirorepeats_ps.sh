#!/usr/bin/bash

set -e

file1=$1

#printf "\tProcessing $file1\n" ;

sed -n '/^\/xstretch/,/^0.1 cm /p' $file1 | 
    grep '^[0-9]' | 
    grep -v -w -e 'setlinewidth' | 
    awk -F' ' 'OFS="\t" {print $2,$3,$4,$7,$8,$9}' | 
    sed -e 's/(//g' -e 's/)//g' 

