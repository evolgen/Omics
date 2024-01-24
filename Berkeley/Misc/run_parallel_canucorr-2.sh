#!/usr/bin/sh

set -e

input=$@

workdir=$(echo $input | awk '{print $1}')
count=$(echo $input | awk '{print $2}')

cd $workdir

number=$(printf "%06d\n" "$count")
if [ ! -e "correction/2-correction/correctReads.${number}.out" ]; then
	cd correction/2-correction/
	echo "Running correction $count"
	./correctReads.sh $count > ./correctReads.${number}.out 2>&1
fi

