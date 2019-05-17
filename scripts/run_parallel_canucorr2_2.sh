#!/usr/bin/sh

set -e

count=$@

cd $workdir

if [ ! -e "correction/2-correction/correctReads.000${count}.out" ]; then
	cd correction/2-correction/
	echo "Running correction $count"
	./correctReads.sh $count > ./correctReads.000${count}.out 2>&1
fi

