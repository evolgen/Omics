#!/usr/bin/bash

# find /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/nucleargenome/SCOs_76/sequence_files/ -type f -name "trim.NT.nex" | parallel -j 32 sh ~/RGP/scripts/selection/parallel_nuclear_hyphy_gard.sh
# ls -d $PWD/sequence_files/*/trim.NT.nex | parallel -j 32 sh ~/RGP/scripts/selection/parallel_nuclear_hyphy_gard.sh

aln1=$@

work_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/nucleargenome/SCOs_76"

sample=$(dirname "${aln1}" | sed -e 's/\/trim.NT.nex$//' -e 's/.*\///'); 

    echo $(dirname $aln1); 
    cd $(dirname $aln1); 

if [ ! -e "${sample}.gard.json" ]; then    
    hyphy GARD --alignment ${aln1} --type codon --tree /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick --output ${sample}.gard_codon.out --output-lf ${sample}.gard_codon.json 1>${sample}.gard_codon.log
    hyphy GARD --alignment ${aln1} --type nucleotide --rv GDD --rate-classes 3 --tree /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick --output ${sample}.gard.out --output-lf ${sample}.gard.json 1>${sample}.gard.log
    printf "  Done\t:\t${sample}\n\n"
fi

