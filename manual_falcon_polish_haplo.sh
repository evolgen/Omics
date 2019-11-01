#!/usr/bin/sh

set -e

cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/NEW/
sh ~/RGP/run_VGP_arrow_haplo.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/umbrosus_pacbio.fastq.gz $PWD/umbrosus_alternate.fasta ~/CONFIGS/falcon/input_bam_umbrosus.fofn
cd arrow_haplo && sh ~/RGP/run_VGP_freebayes.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_umbrosus/S-umbrosus_SEB-10/SEB-10_S144222_R1_001.fastq.gz /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_umbrosus/S-umbrosus_EB-10/SEB-10_S144222_R2_001.fastq.gz $PWD/umbrosus_alternate.t1.fasta

printf "\n\t\tStart - ALASCANUS\n\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/NEW/
sh ~/RGP/run_VGP_arrow_haplo.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/alascanus_pacbio.fastq.gz $PWD/alascanus_alternate.fasta /global/home/users/rohitkolora/CONFIGS/falcon/input_bam_alascanus_new.fofn
printf "\n\t FREEBAYESing\n\n"
cd arrow_haplo && sh ~/RGP/run_VGP_freebayes.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastolobus_alascanus/B-alascanus_SEB-31/Seb31_S134S57_R1_001.fastq.gz /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastolobus_alascanus/B-alascanus_SEB-31/Seb31_S134S57_R2_001.fastq.gz $PWD/alascanus_alternate.t1.fasta

printf "\n\t\tStart - ALEUTIANUS\n\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/NEW/
sh ~/RGP/run_VGP_arrow_haplo.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/aleutianus_pacbio.fastq.gz $PWD/aleutianus_alternate.fasta /global/home/users/rohitkolora/CONFIGS/falcon/input_bam_aleutianus.fofn
printf "\n\t FREEBAYESing\n\n"
cd arrow_haplo && sh ~/RGP/run_VGP_freebayes.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_aleutianus/S-aleutianus_SEB-111/Seb111_S129219_R1_001.fastq.gz /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_aleutianus/S-aleutianus_SEB-111/Seb111_S129219_R2_001.fastq.gz $PWD/aleutianus_alternate.t1.fasta

printf "\n\t\tStart - ENTOMELAS\n\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/purging/NEW/
sh ~/RGP/run_VGP_arrow_haplo.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/entomelas_pacbio.fastq.gz $PWD/entomelas_alternate.fasta /global/home/users/rohitkolora/CONFIGS/falcon/input_bam_entomelas.fofn
printf "\n\t FREEBAYESing\n\n"
cd arrow_haplo && sh ~/RGP/run_VGP_freebayes.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_entomelas/S-entomelas_SEB-8/Seb8_S126220_R1_001.fastq.gz /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_entomelas/S-entomelas_SEB-8/Seb8_S126220_R2_001.fastq.gz $PWD/entomelas_alternate.t1.fasta

printf "\n\t\tStart - PINNIGER\n\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging/NEW/
sh ~/RGP/run_VGP_arrow_haplo.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/pinniger_pacbio.fastq.gz $PWD/pinniger_alternate.fasta /global/home/users/rohitkolora/CONFIGS/falcon/input_bam_pinniger.fofn
printf "\n\t FREEBAYESing\n\n"
cd arrow_haplo && sh ~/RGP/run_VGP_freebayes.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_pinniger/S-pinniger_SEB-72/Seb72_S148217_R1_001.fastq.gz /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_pinniger/S-pinniger_SEB-72/Seb72_S148217_R2_001.fastq.gz $PWD/pinniger_alternate.t1.fasta

printf "\n\t\tStart - ROSACEUS\n\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_rosaceus/06.12.2019/Assembly/purging/NEW/
sh ~/RGP/run_VGP_arrow_haplo.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_rosaceus/06.12.2019/rosaceus_pacbio.fastq.gz $PWD/rosaceus_alternate.fasta /global/home/users/rohitkolora/CONFIGS/falcon/input_bam_rosaceus.fofn 
printf "\n\t FREEBAYESing\n\n"
cd arrow_haplo && sh ~/RGP/run_VGP_freebayes.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_rosaceus/S-rosaceus_SEB-73/Seb-73_S137218_R1_001.fastq.gz /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_rosaceus/S-rosaceus_SEB-73/Seb-73_S137218_R2_001.fastq.gz $PWD/rosaceus_alternate.t1.fasta

printf "\n\t\tStart - RUBERRIMUS\n\n"
cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging/NEW/
sh ~/RGP/run_VGP_arrow_haplo.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/ruberrimus_pacbio.fastq.gz $PWD/ruberrimus_alternate.fasta /global/home/users/rohitkolora/CONFIGS/falcon/input_bam_ruberrimus.fofn
printf "\n\t FREEBAYESing\n\n"
cd arrow_haplo && sh ~/RGP/run_VGP_freebayes.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_ruberrimus/S-ruberrimus_SEB-74/Seb74_S138221_R1_001.fastq.gz /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/Sebastes_ruberrimus/S-ruberrimus_SEB-74/Seb74_S138221_R2_001.fastq.gz $PWD/ruberrimus_alternate.t1.fasta


