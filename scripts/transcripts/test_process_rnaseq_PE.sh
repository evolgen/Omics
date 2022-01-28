#!/usr/bin/sh

set -e 

cd /global/scratch/rohitkolora/Rockfish/Public/Transcriptome/sequencing/SRA/

for file1 in ./S_*/*_1.fastq.gz; do

  file2=$(echo $file1 | sed -e 's/_1.fastq.gz$/_2.fastq.gz/')

echo "  fastp --trim_poly_g --detect_adapter_for_pe --overrepresentation_analysis --thread 20 --length_required 35 --correction --overlap_diff_limit 2 --dont_overwrite --trim_poly_x --disable_quality_filtering --n_base_limit 1 --in1 $file1 --in2 $file2 --out1 ${file1}.adptrm --out2 ${file2}.adptrm"

  outdir=${file1%/*}
  fileout1=$(basename "$file1")
  fileout=$(echo $fileout1 | sed "s/_1.fastq.gz/.flash/")
  filehist=$(echo $file1 | sed -e "s/_1.fastq.gz$/.hist/")

echo "  flash -z -m 30 -x 0.015 -M 75 -t 20 ${file1}.adptrm ${file2}.adptrm -d $outdir -o $fileout "

  fileE=$(echo $fileout | sed -e 's/$/.extendedFrags.fastq.gz/')
  fileR1=$(echo $fileout | sed -e 's/$/.notCombined_1.fastq.gz/')
  fileR2=$(echo $fileR1 | sed -e 's/.notCombined_1.fastq$/.notCombined_2.fastq.gz/')
  fileout1=$(basename "$file1")
  file_F=$(echo $fileE | sed -e 's/.extendedFrags.fastq.gz$/.Forward.fastq.gz/')
  file_R=$(echo $fileE | sed -e 's/.extendedFrags.fastq.gz$/.Reverse.fastq.gz/')

echo "  zcat $fileR1 $fileE | gzip > $file_F && rm $fileR1 $fileE"
echo "  mv $fileR2 $file_R"
done

