#!/usr/bin/sh

set -e

source /global/scratch/rohitkolora/miniconda3/etc/profile.d/conda.sh
conda activate purge_haplotigs

#/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/B_alascanus/05.28.2019/Assembly/purging/WTDBG2
#/global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_aleutianus/05.28.2019/Assembly/purging/WTDBG2


cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/purging/WTDBG2
minimap2 -t 16 -ax map-pb /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/polish/wtdbg2_002.V2.ctg.fa.3.fasta /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/entomelas_pacbio.fastq.gz | samtools view -hF 256 - | samtools sort -@ 16 -m 1G -o aligned.bam -T ./tmp.ali && samtools index aligned.bam && purge_haplotigs readhist -b aligned.bam -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/polish/wtdbg2_002.V2.ctg.fa.3.fasta -t 32
purge_haplotigs readhist -b aligned.bam -t 24 -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/polish/wtdbg2_002.V2.ctg.fa.3.fasta
purge_haplotigs contigcov -i aligned.bam.gencov -l 6 -m 42 -h 190 -o coverage_stats.csv -j 80 -s 80
purge_haplotigs purge -c coverage_stats.csv -t 24 -o entomelas_purged -b aligned.bam -d -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_entomelas/04.04.2019/Assembly/polish/wtdbg2_002.V2.ctg.fa.3.fasta


mkdir -p /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/WTDBG2 /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/FALCON
cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/purging/WTDBG2
minimap2 -t 16 -ax map-pb /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/freebayes/umbrosus_wtdbg2_007.V1.ctg.fa.3.frby.fasta /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/umbrosus_pacbio.fastq | samtools view -hF 256 - | samtools sort -@ 8 -m 1G -o aligned.bam -T /clusterfs/genomicdata/purge_tmp.ali
samtools index aligned.bam
purge_haplotigs readhist -b aligned.bam -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/freebayes/umbrosus_wtdbg2_007.V1.ctg.fa.3.frby.fasta -t 24
purge_haplotigs contigcov -i aligned.bam.gencov -l 15 -m 84 -h 190 -o coverage_stats.csv -j 80 -s 80
purge_haplotigs purge -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/polish/freebayes/umbrosus_wtdbg2_007.V1.ctg.fa.3.frby.fasta -c coverage_stats.csv -t 24 -o umbrosus_purged -b aligned.bam -d


mkdir -p /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging/WTDBG2 /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging/FALCON
cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/purging/WTDBG2
minimap2 -t 16 -ax map-pb /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/polish/freebayes/pinniger_wtdbg2_001.V1.ctg.fa.5.frby.fasta /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/pinniger_pacbio.fastq | samtools view -hF 256 - | samtools sort -@ 8 -m 1G -o aligned.bam -T /clusterfs/genomicdata/purge_tmp.ali
samtools index aligned.bam
purge_haplotigs readhist -b aligned.bam -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/polish/freebayes/pinniger_wtdbg2_001.V1.ctg.fa.5.frby.fasta -t 24
purge_haplotigs contigcov -i aligned.bam.gencov -l 15 -m 84 -h 190 -o coverage_stats.csv -j 80 -s 80
purge_haplotigs purge -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/polish/freebayes/pinniger_wtdbg2_001.V1.ctg.fa.5.frby.fasta -c coverage_stats.csv -t 24 -o pinniger_purged -b aligned.bam -d


mkdir -p /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging/WTDBG2 /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging/FALCON
cd /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/purging/WTDBG2
minimap2 -t 16 -ax map-pb /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/polish/freebayes/ruberrimus_wtdbg2_001.V2.ctg.fa.4.frby.fasta /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/ruberrimus_pacbio.fastq | samtools view -hF 256 - | samtools sort -@ 8 -m 1G -o aligned.bam -T /clusterfs/genomicdata/purge_tmp.ali
samtools index aligned.bam
purge_haplotigs readhist -b aligned.bam -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/polish/freebayes/ruberrimus_wtdbg2_001.V2.ctg.fa.4.frby.fasta -t 24
purge_haplotigs contigcov -i aligned.bam.gencov -l 17 -m 86 -h 190 -o coverage_stats.csv -j 80 -s 80
purge_haplotigs purge -g /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_ruberrimus/04.29.2019/Assembly/polish/freebayes/ruberrimus_wtdbg2_001.V2.ctg.fa.4.frby.fasta -c coverage_stats.csv -t 24 -o ruberrimus_purged -b aligned.bam -d

