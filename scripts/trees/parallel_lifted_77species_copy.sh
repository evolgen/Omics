#!/usr/bin/bash

#set -e -o pipefail

cds=$@

module load samtools

if [ "$#" -ne 1 ]; then
    echo "Need atleast one argument"
    exit 1;
fi

#pep=$(echo $cds | sed -e 's/seq.FNA/seq.FAA/'); 
name=$(dirname "$cds" | sed -e 's/.*\///'); 
new_location=$(dirname "$cds" | sed -e 's/\/Lifted\//\/Lifted\/Sebastes_77\//');

printf $name"\t"$new_location"\n" ; 

if [[ ! -f "${new_location}/seq.FNA.fai" ]]; then
    mkdir -p $new_location ;

    sed -e '/^>/ s/.FUN.*//' $cds >${cds}.edit && samtools faidx ${cds}.edit ;
    cat ~/RGP/scripts/selection/sebastes_species.alascanus.list |
        xargs samtools faidx ${cds}.edit >${new_location}/seq.FNA && samtools faidx ${new_location}/seq.FNA ; 
    
    if [[ -f "${new_location}/seq.FNA.fai" ]]; then
        lines=$(cat ${new_location}/seq.FNA.fai | wc -l);
        if [[ "$lines" -ne 77 ]]; then 
            mv ${new_location}/seq.FNA ${new_location}/seq.FNA.prob; 
        fi
    fi    
fi 

#sed -e '/^>/ s/.FUN.*//' $pep >${pep}.edit && samtools faidx ${pep}.edit ;
#cat ~/RGP/scripts/selection/sebastes_species.alascanus.list |
#    xargs samtools faidx ${pep}.edit >${new_location}/seq.FAA ;


