#!/usr/bin/sh

set -e

#count=$@

###	for count in `seq 49 99`; do echo $count; done | parallel -j 20 sh run_parallel_canucorr1.sh	###

### sh ~/SCRIPTS/run_parallel_canucorr1.sh /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_pinniger/05.08.2019/Assembly/canu/run1/ 49 99

workdir=$1
start=$2
stop=$3

cd $workdir

if [ ! -e "correction/2-correction/correctReads.0000${count}.out" ]; then
	cd correction/2-correction/
	echo "Running correction $count"
	./correctReads.sh $count > ./correctReads.0000${count}.out 2>&1
fi


