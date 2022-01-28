#!/usr/bin/bash

# ls -d $PWD/sequence_files/*/trim.NT.nex | parallel -j 32 sh run_parallel_dndstree_fromaln.sh

aln1=$@

work_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/nucleargenome/SCOs_76"

sample=$(dirname "${aln1}" | sed -e 's/\/trim.NT.nex$//' -e 's/.*\///'); 

    echo $(dirname $aln1); 
    cd $(dirname $aln1); 
    bash ~/RGP/scripts/selection/hyphy_dnds.run.sh Universal $aln1 /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick $(dirname $aln1 | sed -e "s/$/\/${sample}/") ; 


