#!/usr/bin/bash

set -e

for file1 in /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/*minia*; do

  if [[ "$file1" =~ "nigrocinctus" ]]; then
    continue;
  fi
  tmp_path=$(dirname "$file1")
  species=$(echo $file1 | sed -e 's/.*\/minia_//' -e 's/_k.*//' -e 's/_.*//' -e 's/S-/Sebastes_/' -e 's/B-/Sebastolobus_/' -e 's/H-/Helicolenus_/');
  sp_sample=$(echo $file1 | sed -e 's/.*\/minia_//' -e 's/_k.*//');
#  printf "${file1}\t${species}\t${sp_sample}\t${tmp_path}/output/${species}/${sp_sample}/assembly"; echo;
  mv ${file1} ${tmp_path}/output/${species}/${sp_sample}/assembly;


done
  
