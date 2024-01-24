#!/usr/bin/bash

set -e

###sh ~/SCRIPTS/run_salsa4hic.sh \
###             /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbroosus/01.25.2019/Assembly/polish/hic/arima
###             /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/hiC_Arima/SEB_10/ \
###             GATC,GANTC \                                                                                                                    #DNASE
###             /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/wtdbg2_007.V1.ctg.fa \
###             sebumb1
#####

if [ $# -ne 4 ]; then
    printf "\n\n\tPlease provide  parameters as follows\n";
    printf "\t\trun_salsa4hic.sh WORKDIR HiC_FASTQ_PATH REF_FASTA NAME\n\n";
    exit 1;
fi


presentdir=$PWD
workdir=$1
arima_fastq_path=$2
#enzyme_sites=$3
reference_fasta=$3
naming=$4

source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh;
conda activate py2.7
module load java samtools/1.8 minimap2 bwa bedtools kentutils bowtie2 gcc zlib boost salsa

mkdir -p ${workdir} && cd ${workdir}

printf "\n\nYour parameters are as follows:\n";
printf "\tReference = ${reference_fasta}\n";
printf "\tDNase fastq path = ${arima_fastq_path}\n";
printf "\tRestriction enzyme sites = DNase\n";
printf "\tWorking directory = $workdir\n";
printf "\tGiven name for reference = $naming\n";

if [[ ! -e "./deduplicated/repli.bed" ]]; then
    printf "\tRuning the Mapping step for Arima\n";
    run_mapping_arima.sh ${reference_fasta} ${arima_fastq_path} ${naming}
    printf "\n\tConverting the bam to bed\n";
    bamToBed -i ./deduplicated/repli.bam >./deduplicated/repli.bed
    printf "\tSorting the bed\n";
    sort -k 4 -T /clusterfs/genomicdata/ ./deduplicated/repli.bed > /clusterfs/genomicdata/${naming}.tmp && mv /clusterfs/genomicdata/${naming}.tmp ./deduplicated/repli.bed && rm ./deduplicated/repli.bam
fi    

        
if [ ! -e "${reference_fasta}.fai" ]; then
    samtools faidx ${reference_fasta};
fi

printf "\tScaffolding with SALSA\n";
python /global/home/users/rohitkolora/local_modules_sw/salsa/2.2/run_pipeline.py -a ${reference_fasta} -c 1000 -o scaffolds_nomisasbmcor -l ${reference_fasta}.fai -e DNASE -i 3 -b ./deduplicated/repli.bed -m no 1>>LOG 2>>LOG

printf "\tConverting the hic files\n";
convert.sh scaffolds 1>>LOG 2>>LOG
awk -f /global/home/users/rohitkolora/local_modules_sw/3ddna/utils/generate-assembly-file-from-fasta.awk ./scaffolds/scaffolds_FINAL.fasta >./scaffolds/scaffolds_FINAL.assembly


conda deactivate && conda activate py3;
module load busco/3.1 augustus/2.5.5 hmmer blast/2.2.26 gcc;
AUGUSTUS_CONFIG_PATH="/global/home/users/rohitkolora/local_modules_sw/augustus/2.5.5/config/";
mkdir -p ../BUSCO && cd ../BUSCO;
for file1 in ${workdir}/scaffolds/scaffolds_FINAL.fasta; do
    run_name=$(echo $file1 | sed -e 's/^\/.*\/salsa_/salsa_/' -e 's/\/.*//');
    python /global/home/users/rohitkolora/local_modules_sw/busco/3.1/scripts/run_BUSCO.py \
        -i $file1 -c 32 -o ${run_name}_vert \
        -t ./run_${run_name} -m geno -sp zebrafish \
        -l /global/scratch2/rohitkolora/databases/busco/vertebrata_odb9;
done

