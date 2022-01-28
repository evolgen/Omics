#!/usr/bin/bash

set -e

printf "Group\tType\tK\tP\n" #>filt.relax_all.KP.txt; 

find /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/ -type f -name "relax_long.filt.log" | 
    fgrep -e 'REVISE/RELAX_raw/relax_'| 
    while read file1; do 
        name=$(dirname "${file1}" | sed -e 's/.*\/Ortho/Ortho/' -e 's/\/.*//'); 
        printf $name"\tlong\t"; 
        tac $file1 | fgrep -m 1 -e 'Relaxation/intensification parameter (K)' | sed -e 's/.*(K).* //' | tr "\n" "\t"; 
        tac $file1 | fgrep -m 1 -e 'Likelihood ratio test **p =' | sed -e 's/.* //' -e 's/\*\*\.//' | tr "\n" "\t" | sed -e 's/\t$//'; 
        echo; 
    done

find /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/ -type f -name "relax_short.filt.log" | 
    fgrep -e 'REVISE/RELAX_raw/relax_'| 
    while read file1; do 
        name=$(dirname "${file1}" | sed -e 's/.*\/Ortho/Ortho/' -e 's/\/.*//'); 
        printf $name"\tshort\t"; 
        tac $file1 | fgrep -m 1 -e 'Relaxation/intensification parameter (K)' | sed -e 's/.*(K).* //' | tr "\n" "\t"; 
        tac $file1 | fgrep -m 1 -e 'Likelihood ratio test **p =' | sed -e 's/.* //' -e 's/\*\*\.//' | tr "\n" "\t" | sed -e 's/\t$//'; 
        echo; 
    done 

