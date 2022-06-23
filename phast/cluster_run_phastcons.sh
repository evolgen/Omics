#!/bin/bash

#$ -N phastcons
#$ -S /bin/bash
#$ -l h_rt=200:00:00
#$ -l h_vmem=18G
#$ -l highmem

#$ -pe smp 18 

#$ -o /work/kolora/Genome/Multiz/log_$JOB_NAME-$JOB_ID
#$ -j y

export PATH=$PATH:/home/kolora/phast-1.3/bin
cd /work/kolora/Genome/Multiz 

# Processing chunk $ SGE_TASK_ID ...


ls new_mod_mafs/*.maf | parallel -j 18 sh run_new_phastcons_1.sh
ls new_Chunks/*.ss | parallel -j 18 sh run_new_phastcons_2.sh
ls new_Trees/*.cons.mod > Lacerta_cons.txt
ls new_Trees/*.ls new_Chunks/*.ss | parallel -j 18 sh run_new_phastcons_4.sh


