#!/usr/bin/bash

set -e -o pipefail

if [ "$#" -ne 3 ] ; then
    echo "Please run as : sh extract_nonoverlapping_sizes.sh <file> <runname> <size>";
    echo "<size> can be more than 500";
    exit 1;
fi

file=$1
runname=$2
sizes=$3

if [ "${sizes}" -lt 500 ] ; then
    echo "<size> should be a number more than 500";
    exit 1;
fi    

grep -e '^>' $file | sed -e 's/^>//' | grep -v -e 'Sebastes_nigrocinctus' >species_names.txt

module load samtools gcc java raxml
samtools faidx $file
max_size=$(head -n 1 ${file}.fai | cut -f2);
#total_iter=$((max_size-1))
steps=$((sizes-1))
total_iter=$((max_size/sizes))
echo "Total iterations = $total_iter"; #sleep 5

mkdir -p ./${runname}_${sizes};
count=0;
for size in $(seq 1 ${sizes} ${max_size}); do
    count=$((count+1));
    if  [[ "${count}" -lt "${total_iter}" ]]; then
        printf "" >./${runname}_${sizes}/seq_${count}.pranktrim.FAS;
        cat species_names.txt |
            while read speciesname; do
                starts=$size;
                ends=$((starts+steps));
                echo "$file ${speciesname}:${starts}-${ends}"
                samtools faidx $file ${speciesname}:${starts}-${ends} | 
                    sed -e '/^>/ s/:.*//' >>./${runname}_${sizes}/seq_${count}.pranktrim.FAS;
            done
        echo "Done with $count";
        current_size=$((ends-starts));
        cd ./${runname}_${sizes}/;
        raxmlHPC-SSE3 -T 1 -f a -m GTRGAMMA -p $count -x $count -# 10 -s seq_${count}.pranktrim.FAS -n ${runname}_${count} -o Sebastolobus_alascanus &
        cd ../;
    fi
done

#cat ./${runname}_${sizes}/*Branch* >${runname}_${sizes}_raxml_branchlengths.txt
