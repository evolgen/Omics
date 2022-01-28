#!/bin/bash

set -e

module load gcc openmpi transdecoder trinotate

work_dir="/global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA"
cd ${work_dir}

for file in ${work_dir}/star/S_*/*/trinity_out_dir/Trinity-GG.fasta; do

   cd $(dirname "${file}")
   printf "Running orf-finding\n";
   TransDecoder.LongOrfs -m 30 -t ${file/%.gz/}
   printf "Running predictions\n";
   TransDecoder.Predict -t ${file/%.gz/}

   cd ${work_dir}

done


