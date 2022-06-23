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



ls new_Chunks/Lvir_*.ss | parallel -j 18 sh run_new_phylop_1.sh
cat new_phylop_elements/Lvir_*.element-scores.txt | fgrep -v '#chr' | sort -k1,1V -k2,2n >Lacerta_phylop_element-scores.txt
sed -i -e 's/lacvir1\.//' -e '1s/^/chr\tstart\tend\tname\tderiv\tteststat\tpval\n/' Lacerta_phylop_element-scores.txt
ls new_mod_mafs/Lvir_*.maf | parallel -j 18 sh run_new_phylop_2.sh
cat new_phylop_scores/*.scores.wig >Lacerta_phylop_scores.wig
ls new_Chunks/Lvir_*.ss | parallel -j 18 sh run_new_phylop_1_byb.sh
cat new_phylop_elements_bb/*.txt | fgrep -v '#chr' | sort -k1,1V -k2,2n >Lacerta_phylop_base-scores.wig



