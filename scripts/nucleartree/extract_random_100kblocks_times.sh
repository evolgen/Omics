#!/usr/bin/bash

set -e -o pipefail

if [ "$#" -ne 3 ] ; then
    echo "Please run as : sh extract_random_sizes.sh <file> <runname> <size>";
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

rm -f ./random_${runname}_${sizes}/LOG.log;

module load samtools gcc java raxml
samtools faidx $file
max_size=$(head -n 1 ${file}.fai | cut -f2);
limit=$((max_size-sizes))
steps=$((sizes-1))
total_iter=$((max_size/sizes))

mkdir -p ./random_${runname}_${sizes};
echo "" >./random_${runname}_${sizes}/LOG.log;

count=0;
for size in $(seq 0 1 100); do
    count=$((count+1));
    printf "Starting block \t ${count}\n";
        mkdir -p ./random_${runname}_${sizes}/create_${count}/;
        for run in $(seq 1 1 10); do
            printf " $run ";
            starts=$(shuf -i 1-${limit} -n 1);
            ends=$((starts+steps));
            printf "${size}\t${starts}\t${ends}\n" >>./random_${runname}_${sizes}/LOG.log;
            cat species_names.txt |
                while read speciesname; do
                    samtools faidx $file ${speciesname}:${starts}-${ends} | 
                        grep -v '^>' | tr "\n" " " | 
                        sed -e 's/ //g' >>./random_${runname}_${sizes}/create_${count}/${speciesname}.fasta;
               done
            done
            echo;
done

count=0;
for size in $(seq 0 1 100); do
    count=$((count+1));
    printf "" >./random_${runname}_${sizes}/rand_seq_${count}.pranktrim.FAS;
    for files in ./random_${runname}_${sizes}/create_${count}/*.fasta; do
        species_name=$(echo $files | sed -e 's/.*\///' -e 's/\.fasta//');
        printf ">${species_name}\n" >>./random_${runname}_${sizes}/rand_seq_${count}.pranktrim.FAS;
        cat ${files} >>./random_${runname}_${sizes}/rand_seq_${count}.pranktrim.FAS;
        printf "\n" >>./random_${runname}_${sizes}/rand_seq_${count}.pranktrim.FAS;
    done    
    cd ./random_${runname}_${sizes}/;
    raxmlHPC-SSE3 -T 1 -f a -m GTRGAMMA -p $count -x $count -# 10 -s rand_seq_${count}.pranktrim.FAS -n random_${runname}_${count} -o Sebastolobus_alascanus &
    cd ../;
done

