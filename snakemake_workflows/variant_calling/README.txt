
ln002.brc: $  ls -tr /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/fastq/*_*/*/*_R1_* | sed -e 's/\/global.*\/fastq\//"/' -e 's/_R1_.*gz/",/' >>/global/home/users/rohitkolora/RGP/snakemake_workflows/variant_calling/pathchange_config.2nd.json


