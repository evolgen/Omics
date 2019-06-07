#!/usr/bin/sh

set -e

reference_genome="/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/freebayes/umbrosus_wtdbg2_007.V1.ctg.fa.3.frby.fasta"
workdir="/global/scratch/rohitkolora/Rockfish/Genomes/alignments/variantcalling/"

for file1 in /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_[a-m]*/*/*_R1_001.fastq.gz; do

    file2=$(echo $file1 | sed -e 's/_R1_001.fastq.gz/_R2_001.fastq.gz/');
    identity=$(echo $file1 | sed -e 's/.*\///' -e 's/_R1_001.fastq.gz//' -e 's/_.*//');
    species=$(echo $file1 | sed -e "s/\/global\/scratch\/rohitkolora\/Rockfish\/Genomes\/alignments\/variantcalling\///" -e "s/${identity}/ /" -e 's/ .*//');
    genus=$(echo $species | sed -e 's/Sebastes_/S-/' -e 's/Sebastolobus_/B-/');
    samples="${genus}_${identity}";
    
    cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq;
    mkdir -p calls && mkdir -p calls/${species}; 
    mkdir -p calls/${species}/${genus}_${identity};
    mkdir -p calls/${species}/${genus}_${identity}/${identity};
    cd calls/${species}/${genus}_${identity}/${identity};
    module load minimap2 samtools bcftools;
    minimap2 -t 10 -ax sr ${reference_genome} \
        ${file1} ${file2} | \
        samtools view -h - | \
        samtools sort -@ 20 -m 1G \
            -T /clusterfs/genomicdata/map_${identity} \
            -o sort.bam;
    samtools index sort.bam;
    source /global/scratch/rohitkolora/miniconda3/etc/profile.d/conda.sh;
    conda activate gatk;
    /global/scratch/rohitkolora/Software/gatk-4.1.2/gatk \
        --java-options "-Xmx350G" MarkDuplicates \
        --TMP_DIR /clusterfs/genomicdata/ \
        --VALIDATION_STRINGENCY LENIENT \
        --INPUT sort.bam \
        --OUTPUT dedup.bam \
        --METRICS_FILE sort.bam.metrics \
        --MAX_FILE_HANDLES 15000;

done
        
