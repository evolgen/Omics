#!/bin/bash
#SBATCH --job-name=ragoo2umbrosus
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
#SBATCH --output=job_ragoo_%j.log
#SBATCH --error=job_ragoo_%j.err
#
## Command(s) to run (example):

set -e

source /global/scratch/rohitkolora/miniconda3/etc/profile.d/conda.sh
module load minimap2 samtools gcc ragoo survivor

working_dir="/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/comparative/assemblies"
ref_fasta="/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/polish/wtdbg2_003.ctg.V3.fa"
cd $working_dir

for file1 in /global/scratch/rohitkolora/Rockfish/Public/Genome/assemblies/S_*/*.fasta.gz; do 
  
  ref_name=$(echo $ref_fasta | sed -e 's/.*\///')
  species=$(echo $file1 | sed -e 's/.*\/S_/S_/' -e 's/\/.*//')
  name=$(echo $file1 | sed -e 's/.*\///' -e 's/\.gz$//')
  mkdir -p $species && cd $species

  touch ${ref_name}
  if  [[ ${ref_fasta} == *.gz ]] ; then
    zcat $file1 >${ref_name}
  else
     ln -sf ${ref_fasta} ${working_dir}/${species}/${ref_name}
  fi 
  
  touch $name
  if  [[ $file1 == *.gz ]] ; then
    zcat $file1 >$name
  else
     ln -sf $file1 ${working_dir}/${species}/$name
  fi
  
  python3 /global/scratch/rohitkolora/miniconda3/bin/ragoo.py -t 20 -s -m /global/scratch/rohitkolora/miniconda3/bin/minimap2 $name ${ref_name}
  for file1 in ./ragoo_output/pm_alignments/assemblytics_out.Assemblytics_structural_variants.bed; do awk -F'\t' '$NF>=0.75 || $(NF-1)>=0.75' $file1 | cut -f1-11 >${file1}.11col; SURVIVOR convertAssemblytics ${file1}.11col 30 ${file1/%.bed/.vcf}; done
  
  if  [[ ${ref_fasta} == *.gz ]] ; then
    rm ${ref_name}
  fi
  if  [[ ${file1} == *.gz ]] ; then
    rm $name
  fi

  cd $working_dir

done


