#! /bin/bash

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

sh run_mapping_arima.sh ${reference_fasta} ${arima_fastq_path} ${naming} 2>>${naming}/LOG

IN_DIR=$2
WORK_DIR="$PWD/"$3
FASTQ_R1=$(ls $IND_DIR | grep '_R1_')
FASTQ_R2=$(ls $IND_DIR | grep '_R2_')

SRA=$(ls $IND_DIR | grep '_R1_' | sed -e 's/.*\///' -e 's/_R1_.*//')
LABEL="salsa"

BWA='/global/home/groups/consultsw/sl-7.x86_64/modules/bwa/0.7.17/bwa'
MINIMAP2='/global/home/users/rohitkolora/local_modules_sw/minimap2/2.15/minimap2'
SAMTOOLS='/global/home/groups/consultsw/sl-7.x86_64/modules/samtools/1.8/bin/samtools'
REF=$1
FAIDX=$1".fai"
PREFIX=$3
RAW_DIR=$WORK_DIR"/raw"
FILT_DIR=$WORK_DIR"/filtered"
FILTER='/global/home/users/rohitkolora/local_modules_sw/salsa/2.2/filter_five_end.pl'
COMBINER='/global/home/users/rohitkolora/local_modules_sw/salsa/2.2/two_read_bam_combiner.pl'
STATS='/global/home/users/rohitkolora/local_modules_sw/salsa/2.2/get_stats.pl'
PICARD='/global/home/users/rohitkolora/local_modules_sw/picard/2.19.0/picard.jar'
TMP_DIR=$WORK_DIR"/temp"
PAIR_DIR=$WORK_DIR"/paired"
REP_DIR=$WORK_DIR"/deduplicated"
REP_LABEL=$LABEL\_rep1
MERGE_DIR=$WORK_DIR"/merged"
MAPQ_FILTER=20
CPU=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
MEMORY=$(free -g | tail -n +2 | head -n 1 | awk '{print $4}' | sed -e 's/^/-Xmx/' -e 's/$/G/')

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $WORK_DIR ] || mkdir -p $WORK_DIR
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
###$BWA index -a bwtsw -p $PREFIX $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
###$BWA mem -t $CPU $REF $IN_DIR/$SRA\_1.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_1.bam
$MINIMAP2 -t $CPU $REF $FASTQ_R1 | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
###$BWA mem -t $CPU $REF $IN_DIR/$SRA\_2.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_2.bam
$MINIMAP2 -t $CPU $REF $FASTQ_R2 | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/$SRA\_1.bam $FILT_DIR/$SRA\_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA.bam -

echo "### Step 3.B: Add read group"
java $MEMORY -Djava.io.tmpdir=$TMP/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#	REP_NUM=X #number of the technical replicate set e.g. 1
#	REP_LABEL=$LABEL\_rep$REP_NUM
#	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#	java -Xmx300G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java $MEMORY -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
#
#	INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as biological replicates
#   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
#
#	java -Xmx300G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
#
#	$SAMTOOLS index $MERGE_DIR/$LABEL.bam

# perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# echo "Finished Mapping Pipeline through merging Biological Replicates"
