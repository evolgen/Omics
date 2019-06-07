#!/bin/bash
#SBATCH --job-name=CLAIRVOYANTE
#
#SBATCH --account=fc_genomicdata
#
# QoS: must be savio_long for jobs > 3 days
#
# Partition:
#SBATCH --partition=savio
#
# Number of tasks needed for use case (example):
#SBATCH -N 1
#SBATCH -c 20
#
# Wall clock limit (7 days in this case):
#SBATCH --time=72:00:00
#SBATCH --output=/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/Variants/job_clairvoyante_%j.log
#SBATCH --error=/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/Variants/job_clairvoyante_%j.err

#
## Command(s) to run (example):

module load samtools java vcftools
cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/Variants/clairvoyante; 
mkdir -p clairvoyante; 

source /global/scratch/rohitkolora/miniconda3/etc/profile.d/conda.sh
conda activate clairvoyante

# samtools faidx /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/polish/wtdbg2_003.ctg.V3.fa
samtools index /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/Variants/clairvoyante/wtdbg2_003.ctg.V3.srt.bam
awk '{print $1}' /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/polish/wtdbg2_003.ctg.V3.fa.fai | while read ctg; do 
    python /global/home/users/rohitkolora/local_modules_sw/Clairvoyante/clairvoyante/callVarBam.py --threads 20 --bam_fn /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/Variants/clairvoyante/wtdbg2_003.ctg.V3.srt.bam --ref_fn /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/wtdbg2/polish/wtdbg2_003.ctg.V3.fa --minCoverage 8 --threshold 0.2 --call_fn ./clairvoyante/${ctg}.vcf --chkpnt_fn /global/home/users/rohitkolora/local_modules_sw/Clairvoyante/trainedModels/fullv3-pacbio-ngmlr-hg001-hg19/learningRate1e-3.epoch999 --ctgName $ctg; 
done

grep '^#' clairvoyante/ctg1.vcf >wtdbg2_003.V3.clairvoyante_selfmap_pacbio.vcf; cat clairvoyante/*.vcf | grep -v '^#' >tmp; 
sort -k1,1V -k2,2n tmp >>wtdbg2_003.V3.clairvoyante_selfmap_pacbio.vcf; rm tmp




