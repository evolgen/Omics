#!/bin/bash

#$ -N phylop
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

ls new_Chunks/Lvir_*.ss | parallel -j 18 sh run_new_phylop_2_byb.sh
cat new_phylop_elements_bb/*.txt | fgrep -v '#chr' | sort -k1,1V -k2,2n >Lacerta_phylop_base-scores.txt



