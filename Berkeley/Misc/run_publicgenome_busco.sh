#!/usr/bin/bash

set -e

module load busco/3.1 augustus/2.5.5 hmmer blast/2.2.26 gcc
AUGUSTUS_CONFIG_PATH="/global/home/users/rohitkolora/local_modules_sw/augustus/2.5.5/config/"

for file1 in /global/scratch2/rohitkolora/Rockfish/Public/Genome/assemblies/S_*/*.fasta.gz; do

#	zcat $file1 >${file1/%.gz/};
	direcname=$(dirname $file1);
	name1=$(dirname $file1 | sed -e 's/\/$//' -e 's/.*\///');
	cd $direcname;
#	python /global/home/users/rohitkolora/local_modules_sw/busco/3.1//scripts/run_BUSCO.py -i ${file1/%.gz/} -c 24 -o ${name1}_vert -m geno -sp zebrafish -l /global/scratch2/rohitkolora/databases/busco/vertebrata_odb9
	python /global/home/users/rohitkolora/local_modules_sw/busco/3.1//scripts/run_BUSCO.py -i ${file1/%.gz/} -c 24 -o ${name1}_acti -m geno -sp zebrafish -l /global/scratch2/rohitkolora/databases/busco/actinopterygii_odb9
	python /global/home/users/rohitkolora/local_modules_sw/busco/3.1//scripts/run_BUSCO.py -i ${file1/%.gz/} -c 24 -o ${name1}_cvg -m geno -sp zebrafish -l /global/scratch2/rohitkolora/databases/busco/CVG

done

