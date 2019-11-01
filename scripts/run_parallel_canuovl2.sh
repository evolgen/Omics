#!/usr/bin/sh

set -e

count=$@

cd /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/pacbio/S_umbrosus/01.25.2019/Assembly/canu/run1/trimming/1-overlapper/

if [ ! -e "overlap.000${count}.out" ]; then
	echo "Running overlap $count"
	./overlap.sh $count > ./overlap.000${count}.out 2>&1
fi

