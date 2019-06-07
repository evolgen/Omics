#!/usr/bin/sh

set -e

###sh ~/SCRIPTS/run_salsa4hic.sh \
###             /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/hic/arima
###             /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/hiC_Arima/SEB_10/ \
###             GATC,GANTC \                                                                                                                    #DNASE
###             /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/wtdbg2_007.V1.ctg.fa \
###             sebumb1
#####

if [ $# -ne 4 ]; then
    printf "\n\n\tPlease provide  parameters as follows\n";
    printf "\t\trun_salsa4hic.sh WORKDIR HiC_FASTQ_PATH REF_FASTA NAME\n\n";
    exit 1;
fi


workdir=$1
arima_fastq_path=$2
#enzyme_sites=$3
reference_fasta=$3
naming=$4

source /global/scratch/rohitkolora/miniconda3/etc/profile.d/conda.sh;
conda activate py2.7
module load java samtools/1.8 minimap2 bwa bedtools kentutils bowtie2 gcc zlib boost salsa

mkdir -p ${workdir} && cd ${workdir}
printf "\n\nYour parameters are as follows:\n";
printf "\tReference = ${reference_fasta}\n";
printf "\tArima fastq path = ${arima_fastq_path}\n";
#printf "\tRestriction enzyme sites = $enzyme\n";
printf "\tWorking directory = $workdir\n";
printf "\tGiven name for reference = $naming\n";

printf "\tRuning the Mapping step for Arima\n\t";
echo "run_mapping_arima.sh ${reference_fasta} ${arima_fastq_path} ${naming}"
sh run_mapping_arima.sh ${reference_fasta} ${arima_fastq_path} ${naming}

printf "\tConverting the bam to bed\n";
bamToBed -i ${naming}/deduplicated/repli.bam >${naming}/deduplicated/repli.bed 1>${naming}/LOG 2>>${naming}/LOG
sort -k 4 -T /clusterfs/genomicdata/ ${naming}/deduplicated/repli.bed > /clusterfs/genomicdata/${naming}.tmp && mv /clusterfs/genomicdata/${naming}.tmp ${naming}/deduplicated/repli.bed

if [ ! -e "${reference_fasta}.fai" ]; then
    samtools faidx ${reference_fasta};
fi

printf "\tScaffolding with SALSA\n";
python /global/home/users/rohitkolora/local_modules_sw/salsa/2.2/run_pipeline.py -a ${reference_fasta} -c 1000 -o ${naming}/scaffolding -l ${reference_fasta}.fai -e GATC,GAATC,GAGTC,GACTC,GATTC -i 3 -b ${naming}/deduplicated/repli.bed -m yes 1>>${naming}/LOG 2>>${naming}/LOG

printf "\tConverting the hic files\n";
convert.sh ${naming}/scaffolding 1>>${naming}/LOG 2>>${naming}/LOG


