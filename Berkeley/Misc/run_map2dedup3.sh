#!/usr/bin/sh

set -e

reference_genome="/global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/freebayes/umbrosus_wtdbg2_007.V1.ctg.fa.3.frby.fasta"
workdir="/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/variantcalling/"

#for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_[a-m]*/*/*_R1_001.fastq.gz; do
#for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_[b-mo-z]*/*/*_R1_001.fastq.gz; do
for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebast*/*/*_R1_001.fastq.gz; do

    file2=$(echo $file1 | sed -e 's/_R1_001.fastq.gz/_R2_001.fastq.gz/');
    identity=$(echo $file1 | sed -e 's/.*\///' -e 's/_R1_001.fastq.gz//');# -e 's/_.*//');
    echo $identity
    species=$(echo $file1 | sed -e 's/.*\/illumina\/fastq\///' -e 's/\/.*//');
    echo $species
    genus=$(echo $file1 | sed -e 's/.*\/Sebastes_/S-/' -e 's/.*\/Sebastolobus_/B-/' -e 's/\/.*//');# -e 's/\// /' -e 's/ .*//');
    echo $genus
    samples="${genus}_${identity}";
    echo $samples
    
    cd ${workdir};
    mkdir -p calls && mkdir -p calls/${species}; 
    mkdir -p calls/${species}/${genus}_${identity};
    mkdir -p calls/${species}/${genus}_${identity}/${identity};
    cd calls/${species}/${genus}_${identity}/${identity};
    echo $PWD;
    if [ ! -e "sort.bam.bai" ]; then
        printf "\tMapping for ${species} - ${identity}\n";
        module load minimap2 samtools bcftools;
        minimap2 -t 32 -ax sr ${reference_genome} \
            ${file1} ${file2} | \
            samtools view -h - >/clusterfs/genomicdata/${identity}.bam; #| \
        printf "\tSorting ${species} - ${identity}\n";
        samtools sort -@ 32 -m 1G \
            -T /clusterfs/genomicdata/map_${identity} \
            -o sort.bam /clusterfs/genomicdata/${identity}.bam; 
        rm /clusterfs/genomicdata/${identity}.bam &
        samtools index sort.bam;
    fi
    if [ ! -e "dedup.bam" ]; then
        source /global/scratch2/rohitkolora/miniconda3/etc/profile.d/conda.sh;
        conda activate gatk;
        printf "\tDeduplicating for ${species} - ${identity}\n";
        /global/scratch2/rohitkolora/Software/gatk-4.1.2/gatk \
            --java-options "-Xmx600G" MarkDuplicates \
            --TMP_DIR /clusterfs/genomicdata/rockfish/tmp \
            --VALIDATION_STRINGENCY LENIENT \
            --INPUT sort.bam \
            --OUTPUT dedup.bam \
            --METRICS_FILE sort.bam.metrics \
            --MAX_FILE_HANDLES 25000;
        samtools index dedup.bam;
    fi

done
        
