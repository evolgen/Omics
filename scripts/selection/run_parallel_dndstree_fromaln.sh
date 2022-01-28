#!/usr/bin/bash

# ls -d $PWD/sequence_files/*/NT.FNA.trim.aln | parallel -j 32 sh run_parallel_dndstree_fromaln.sh

aln1=$@

work_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/nucleargenome/BUSCO"

sample=$(dirname "${aln1}" | sed -e 's/\/NT.FNA.trim.aln$//' -e 's/.*\///'); 

if [ ! -e "${work_dir}/logs/${sample}_NT.rundnds.done" ] ; then 
    echo $(dirname $aln1); 
    cd $(dirname $aln1); 
    bash ~/RGP/scripts/selection/hyphy_dnds.run.sh Universal $aln1 /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick $(dirname $aln1 | sed -e "s/$/\/${sample}/") ; 
    printf "" >${work_dir}/logs/${sample}_NT.rundnds.done; 
fi 


