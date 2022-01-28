#!/usr/bin/bash

set -e

#####
###sh ~/RGP/run_hic_juicer-forplot.sh \
###             /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/hic/	
###             /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/Seb10_S144_R1.fastq \	
###             /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/Seb10_S144_R2.fastq \
###             /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/wtdbg2_007.V1.ctg.fa \
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

mkdir -p $workdir && cd $workdir
mkdir -p fastq
printf "\tCreating reference and fastq files\n";
ln -sf $reference $PWD/${name}_ref.fasta;  

ln -sf $fastq1 $PWD/fastq/hic_fastq_R1.fastq.gz  
ln -sf $fastq2 $PWD/fastq/hic_fastq_R2.fastq.gz 

#source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh
#conda activate py2.7;
printf "\tCreating restriction sites file\n";
if [ ! -e "${name}.chrom.sizes" ]; then
    python /global/home/users/rohitkolora/local_modules_sw/juicer/1.6.2/misc/generate_site_positions.py HindIII ${name} $reference >${name}_HindIII.txt
    awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${name}_HindIII.txt >${name}.chrom.sizes;
fi

if [ ! -e "${reference}.fai" ]; then
    samtools faidx ${reference} ;
fi
if [ ! -e "${reference}.0123" ] ; then
    printf "\tIndexing the reference with bwa\n";
    /global/scratch2/rohitkolora/Software/bwa-mem2/bwa-mem2 index ${reference} ;
fi    

printf "\tRunning JUICER to create HiC files \n";
if [ -e "./aligned" ]; then
    rm -fr ./aligned/
fi

cut -f1,2 ${reference}.fai >${name}.length
if [ ! -e "./aligned/merged_nodups.txt" ]; then
    bash ~/local_modules_sw/juicer/1.6.2/CPU/juicer.bwamem2.sh -z ${reference} -d $workdir -s HindIII -p ${name}.length -t $(eval nproc) -g ${name} -y ${name}_HindIII.txt 1>LOG 2>>LOG; 
    rm -f ${name}_ref.fasta.* aligned/alignable.bam aligned/collisions*.bam aligned/mapq0.bam aligned/unmapped.bam splits/*.sam splits/*.bam ;
    rm -f ./aligned/merged_nodups.txt ; 
fi    


