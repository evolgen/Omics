#!/bin/bash
#SBATCH --job-name=rnaseq_PE
#
#SBATCH --account=fc_genomicdata
#
# QoS: must be savio_long for jobs > 3 days
#
# Partition:
#SBATCH --partition=savio
#
# Number of tasks needed for use case (example):
#SBATCH --ntasks=20
#SBATCH --nodes=1
#
# Wall clock limit (7 days in this case):
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/job_rnaseq_PE.log
#SBATCH --error=/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/job_rnaseq_PE.err
#
## Command(s) to run (example):

set -e

module load gcc openmpi fastp flash

cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

for file1 in ./S*/*_1.fastq.gz; do

  file2=$(echo $file1 | sed -e 's/_1.fastq.gz$/_2.fastq.gz/')
  filename=$(echo $file1 | sed -e 's/_1.fastq.gz//')
  
  fastp --trim_poly_g --detect_adapter_for_pe --overrepresentation_analysis --html ${filename}.html -j ${filename}.json --thread 20 --length_required 35 --correction --overlap_diff_limit 2 --trim_poly_x --disable_quality_filtering --n_base_limit 1 --in1 $file1 --in2 $file2 --out1 ${file1}.adptrm --out2 ${file2}.adptrm

  outdir=${file1%/*}
  fileout1=$(basename "$file1")
  fileout=$(echo $fileout1 | sed "s/_1.fastq.gz/.flash/")
  filehist=$(echo $file1 | sed -e "s/_1.fastq.gz$/.hist/")

  flash -z -m 30 -x 0.015 -M 75 -t 20 ${file1}.adptrm ${file2}.adptrm -d $outdir -o $fileout
  rm ${file1}.adptrm ${file2}.adptrm

  fileE=$(echo $fileout | sed -e 's/$/.extendedFrags.fastq.gz/')
  fileR1=$(echo $fileout | sed -e 's/$/.notCombined_1.fastq.gz/')
  fileR2=$(echo $fileR1 | sed -e 's/.notCombined_1.fastq.gz$/.notCombined_2.fastq.gz/')
  fileout1=$(basename "$file1")
  file_F=$(echo $fileE | sed -e 's/.extendedFrags.fastq.gz$/.Forward.fastq.gz/')
  file_R=$(echo $fileE | sed -e 's/.extendedFrags.fastq.gz$/.Reverse.fastq.gz/')

  zcat ${outdir}/$fileR1 ${outdir}/$fileE | gzip > ${outdir}/$file_F
#  rm ${outdir}/$fileR1 ${outdir}/$fileE
  cp ${outdir}/$fileR2 ${outdir}/$file_R

done



