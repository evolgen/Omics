#!/usr/bin/sh

set -e

#count=$@

###	for count in `seq 100 250`; do echo $count; done | parallel -j 20 sh run_parallel_canucorr1_1.sh	###

###sh ~/SCRIPTS/run_parallel_canucorr-1.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/canu/run1/ 1 209

workdir=$1
start=$2
stop=$3

for count1 in `seq ${start} ${stop}`; do 
	printf "$workdir $count1\n";
done | parallel -j 20 sh ~/SCRIPTS/run_parallel_canucorr-2.sh 

