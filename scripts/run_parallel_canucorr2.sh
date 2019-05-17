#!/usr/bin/sh

set -e

#count=$@

###	for count in `seq 100 250`; do echo $count; done | parallel -j 20 sh run_parallel_canucorr2.sh	###

###sh ~/SCRIPTS/run_parallel_canucorr2.sh /global/scratch/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/canu/run1/ 100 209

workdir=$1
start=$2
stop=$3

cd $workdir

if [ ! -e "correction/2-correction/correctReads.000${count}.out" ]; then
	cd correction/2-correction/
	echo "Running correction $count"
	./correctReads.sh $count > ./correctReads.000${count}.out 2>&1
fi

