#!/usr/bin/bash

set -e

cut -f4 list_samples.txt | sort -u | while read sample; do
#grep schlegelii list_samples.txt | cut -f4 | sort -u | while read sample; do 


  echo "Running for $sample"
  grep -w "$sample" list_samples.txt | sort -k4,4 -k2,2 | while read line1; do 
    run=$(echo $line1 | awk -F' ' '{print $1}'); 
    type=$(echo $line1 | awk -F' ' '{print $2}'); 
    species=$(echo $line1 | awk -F' ' '{print $3}'); 
    mkdir -p Trinity/${species}; 

    if [ "$type" == "PAIRED" ]; then 
       printf "\tStarted paired run for $run\n";
       cat ./*/${run}.flash.notCombined_1.fastq >>./Trinity/${species}/${sample}_1.fastq; printf "\t\tCatted for ./*/${run}.flash.notCombined_1.fastq\n";
       cat ./*/${run}.flash.notCombined_2.fastq >>./Trinity/${species}/${sample}_2.fastq; printf "\t\tCatted for ./*/${run}.flash.notCombined_1.fastq\n";
       cat ./*/${run}.flash.extendedFrags.fastq >>./Trinity/${species}/${sample}_tmp.fastq; printf "\t\tCatted for ./*/${run}.flash.extendedFrags.fastq\n";
    else 
       printf "\tStarted single run for $run\n";
       cat ./*/${run}.fastq.gz.adptrm >>./Trinity/${species}/${sample}_tmp.fastq; printf "\t\tCatted for ./*/${run}.fastq.adptrm\n";
    fi
  done 

   printf "  Combine yourself for $sample\n"
#   cat ./Trinity/${species}/${sample}_tmp.fastq >>./Trinity/${species}/${sample}_1.fastq; printf "\t\tJoined singles for ./Trinity/${species}/${sample}_tmp.fastq\n";
#   rm ./Trinity/${species}/${sample}_tmp.fastq;
   printf "Done with $sample\n"

#   gzip ./Trinity/${species}/${sample}_1.fastq 
#   gzip ./Trinity/${species}/${sample}_2.fastq 

done


