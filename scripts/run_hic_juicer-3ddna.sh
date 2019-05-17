#!/usr/bin/bash

set -e

#####
###sh ~/SCRIPTS/run_hic_juicer-3ddna.sh \
###             /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/hic/	
###             /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/Seb10_S144_R1.fastq \	
###             /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/Seb10_S144_R2.fastq \
###             /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/wtdbg2_007.V1.ctg.fa \
###             sebumb1
#####

if [ $# -ne 5 ]; then
printf "\n\n\tPlease provide 4 parameters as follows\n";
printf "\t\trun_hic_juicer-3ddna.sh WORKDIR HiC_FASTQ1 HIC_FASTQ2 REF_FASTA NAME\n\n";
exit 1;
fi

workdir=$1
fastq1=$2
fastq2=$3
reference=$4
name=$5

module load java samtools/1.8 minimap2 bwa bedtools kentutils bowtie2 gcc zlib boost 3ddna 

cd $workdir
mkdir -p fastq
printf "\tCreating reference and fastq files\n";
##ln -sf $reference $PWD/${name}_ref.fasta;  

##ln -sf $fastq1 $PWD/fastq/hic_fastq_R1.fastq.gz  
##ln -sf $fastq2 $PWD/fastq/hic_fastq_R2.fastq.gz 

source /global/scratch/rohitkolora/miniconda3/etc/profile.d/conda.sh
conda activate py2.7;
printf "\tCreating restriction sites file\n";
##python /global/home/users/rohitkolora/local_modules_sw/juicer/1.6.2/misc/generate_site_positions.py HindIII ${name} ${name}_ref.fasta >${name}_HindIII.txt;
##awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${name}_HindIII.txt >${name}.chrom.sizes;

##samtools faidx ${name}_ref.fasta &
printf "\tIndexing the reference with bwa\n";
##bwa index $PWD/${name}_ref.fasta;

printf "\tRunning JUICER pipeline\n";
##~/local_modules_sw/juicer/1.6.2/CPU/juicer.sh -z ${name}_ref.fasta -d $workdir -s HindIII -p ${name}_ref.fasta.fai -t 32 -g ${name} -y ${name}_HindIII.txt 1>LOG 2>>LOG; 
printf "JUICER done - now for further scaffolding\n";

mkdir -p 3ddna && cd 3ddna;
printf "\tRunning 3D-DNA\n";
run-asm-pipeline.sh -q 20 --editor-coarse-resolution 2500000 --editor-coarse-region 12500000 --editor-fine-resolution 100000 -i 1000 ../${name}_ref.fasta ../aligned/merged_nodups.txt 1>>LOG 2>>LOG;
printf "3D-DNA done - Now for the stats\n";

final_assembly=$(echo "${name}_ref.fasta" | sed -e 's/\.fasta$/.FINAL.fasta/')
assembly-stats $final_assembly


