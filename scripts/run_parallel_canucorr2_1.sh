#!/usr/bin/sh

set -e

#count=$@

###	for count in `seq 100 250`; do echo $count; done | parallel -j 20 sh run_parallel_canucorr2.sh	###

###sh ~/SCRIPTS/run_parallel_canucorr2_1.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/canu/run1/ 100 209

workdir=$1
start=$2
stop=$3

for count in `seq $start $stop`; do echo $count; done | parallel -j 20 sh ~/SCRIPTS/run_parallel_canucorr2_2.sh

