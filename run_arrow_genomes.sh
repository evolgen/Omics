#!/usr/bin/sh

set -e

##printf "\tRunning entomelas\n"
##cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/purging/FALCON
##mkdir -p arrow && cd arrow
##sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/entomelas_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/purging/FALCON/entomelas_purged.fasta ~/CONFIGS/falcon/input_bam_entomelas.fofn

##printf "\tRunning umbrosus\n"
##cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/FALCON
##mkdir -p arrow && cd arrow
##sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/umbrosus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/FALCON/umbrosus_purged.fasta ~/CONFIGS/falcon/input_bam_umbrosus.fofn

##printf "\tRunning pinniger\n"
##cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging/FALCON
##mkdir -p arrow && cd arrow
##sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/pinniger_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging/FALCON/pinniger_purged.fasta ~/CONFIGS/falcon/input_bam_pinniger.fofn

##printf "\tRunning ruberrimus\n"
##cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging/FALCON
##mkdir -p arrow && cd arrow
##sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/ruberrimus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging/FALCON/ruberrimus_purged.fasta ~/CONFIGS/falcon/input_bam_ruberrimus.fofn

printf "\tRunning alascanus\n"
cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON
mkdir -p arrow && cd arrow
if [ -e "/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq" ]; then
    sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON/alascanus_purged.fasta ~/CONFIGS/falcon/input_bam_alascanus.fofn;
else
   zcat /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq.gz >/clusterfs/genomicdata/rockfish/tmp/alascanus_pacbio.fastq;
   sh ~/SCRIPTS/run_VGP_arrow.sh /clusterfs/genomicdata/rockfish/tmp/alascanus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/FALCON/alascanus_purged.fasta ~/CONFIGS/falcon/input_bam_alascanus.fofn;
fi    

printf "\tRunning aleutianus\n"
cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON
mkdir -p arrow && cd arrow
if [ -e "/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/aleutianus_pacbio.fastq" ]; then
    sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/aleutianus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON/aleutianus_purged.fasta ~/CONFIGS/falcon/input_bam_aleutianus.fofn;
else
    zcat /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/aleutianus_pacbio.fastq.gz >/clusterfs/genomicdata/rockfish/tmp/aleutianus_pacbio.fastq;
    sh ~/SCRIPTS/run_VGP_arrow.sh /clusterfs/genomicdata/rockfish/tmp/aleutianus_pacbio.fastq /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/FALCON/aleutianus_purged.fasta ~/CONFIGS/falcon/input_bam_aleutianus.fofn;
fi    


