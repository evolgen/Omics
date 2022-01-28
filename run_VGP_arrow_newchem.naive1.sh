#!/usr/bin/bash

set -e

#####
###sh ~/SCRIPTS/run_VGP_arrow.sh \
###		/path/to/input_pacbio.fastq \
###		/path/to/ref.fa \
###     /path/to/input_bam.fofn
#####

if [ $# -ne 3 ]; then
	printf "\n\n\tPlease provide 3 parameters as follows\n";
	printf "\t\trun_VGP_arrow_newchem.naive.sh FASTQ REF_FASTA BAM_FOFN\n\n";
	exit 1;
fi

fastq=$1
fastaseq=$2
bam_fofn=$3

reference=$(echo $fastaseq | sed -e 's/.*\///' -e 's/\.fasta$//' -e 's/\.fas$//' -e 's/\.fa$//' -e 's/\.FASTA$//' -e 's/\.FAS$//' -e 's/\.FA$//')

#source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh
#conda activate pb-assembly
module load minimap2 samtools #smrtlink

count=0
mkdir -p arrow && cd arrow

printf "\tCreating reference file\n";
if [[ ! -e "${reference}.${count}.fasta" || ! -e "${reference}.${count}.fasta.fai" ]]; then
  ln -sf $fastaseq ${reference}.${count}.fasta
  samtools faidx ${reference}.${count}.fasta
fi

uniqID=$(echo $BASHPID)

for iter in 1; do
    if [ -e "${reference}.t1.fasta" ]; then
        printf "\t${reference}.t1.fasta EXISTS, hence skipping arrow polishing\n";
        exit 0;
    fi    
	if [ ! -e "pbmap_${reference}.${count}.bam" ]; then
		pbmm2 align --preset SUBREAD ${reference}.${count}.fasta $fastq pbmap_${reference}.${count}.bam --sort -j 26 -J 6;
#        mv /clusterfs/genomicdata/rockfish/tmp/pbmap_${reference}.${count}.bam .;
        printf "\t pbmap_${reference}.${count}.bam mapped\n";
    fi
    if [ ! -e "pbmap_${reference}.pb.bam" ]; then
        printf "\t pbmap_${reference}.${count}.bam conversion to arrow compatible\n";
        /global/scratch2/rohitkolora/miniconda3/envs/smrtlink/share/smrtlink-tools-7.0.1.66975-0/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/private/pacbio/pbbam/binwrap/pbbamify --input=pbmap_${reference}.${count}.bam --output=pbmap_${reference}.pb.bam ${reference}.${count}.fasta $bam_fofn
        printf "\t pbmap_${reference}.pb.bam modified\n";
	fi	
	if [ ! -e "pbmap_${reference}.pbsort.bam.bai" ]; then
        printf "\t pbmap_${reference}.pbsort.bam sorting\n";
        samtools sort -T ./${uniqID} -@ 32 -m 2G -o pbmap_${reference}.pbsort.bam pbmap_${reference}.pb.bam
		pbindex pbmap_${reference}.pbsort.bam;
        printf "\t pbmap_${reference}.pbsort.bam INDEXED\n";
	fi	
	printf "\tVariant calling pbmap_${reference}.pbsort.bam\n";
#    module load smrtlink/v600;
#    conda activate smrtlink
#    variantCaller pbmap_${reference}.pbsort.bam --algorithm=arrow -j 32 -r ${reference}.${count}.fasta -o ${reference}.t1.fasta;
#	samtools faidx ${reference}.t1.fasta;
	count=$((count+1));
done

printf "\tNeeds variantcalling pbmap to produce ${reference}.t1.fasta\n";

