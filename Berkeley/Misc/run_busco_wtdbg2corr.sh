#!/usr/bin/bash

set -e

###/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/wtdbg2/ruberrimus_wtdbg2_000.V6.ctg.fa

module load busco/3.1 augustus/2.5.5 hmmer blast/2.2.26 gcc
AUGUSTUS_CONFIG_PATH="/global/home/users/rohitkolora/local_modules_sw/augustus/2.5.5/config/";

file1=$@

##ls -1 /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/*/*/Assembly/wtdbg2/*wtdbg2_000.V*.ctg.fa | parallel -j 8 sh ~/RGP/run_busco_wtdbg2corr.sh 
##for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/*/*/Assembly/wtdbg2/*wtdbg2_000.V*.ctg.fa; do

    directory=$(dirname $file1)
    filename=$(basename $file1)
    run_name=$(echo $filename | sed -e 's/\.ctg.fa$//')
    mkdir -p ${directory}/BUSCO && cd ${directory}/BUSCO;
    python /global/home/users/rohitkolora/local_modules_sw/busco/3.1/scripts/run_BUSCO.py -i $file1 -c 4 -o ${run_name}_vert -t ./run_${run_name} -m geno -sp zebrafish -l /global/scratch2/rohitkolora/databases/busco/vertebrata_odb9; 
    cd ../; 

##done

