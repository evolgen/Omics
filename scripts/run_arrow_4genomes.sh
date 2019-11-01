#!/usr/bin/sh

set -e

##printf "\tRunning entomelas\n"
##cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/purging
##mkdir -p arrow && cd arrow
##sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/entomelas_pacbio.fastq /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/purging/entomelas_purged.fasta ~/CONFIGS/falcon/input_bam_entomelas.fofn

printf "\tRunning umbrosus\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging
mkdir -p arrow && cd arrow
sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/umbrosus_pacbio.fastq /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/umbrosus_purged.fasta ~/CONFIGS/falcon/input_bam_umbrosus.fofn

printf "\tRunning pinniger\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging
mkdir -p arrow && cd arrow
sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/pinniger_pacbio.fastq /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging/pinniger_purged.fasta ~/CONFIGS/falcon/input_bam_pinniger.fofn

printf "\tRunning ruberrimus\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging
mkdir -p arrow && cd arrow
sh ~/SCRIPTS/run_VGP_arrow.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/ruberrimus_pacbio.fastq /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging/ruberrimus_purged.fasta ~/CONFIGS/falcon/input_bam_ruberrimus.fofn

