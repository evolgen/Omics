#!/usr/bin/sh

set -e -o pipefail

#ls /global/home/users/rohitkolora/RGP/scripts/orthology/illumina_samples_refgen.txt
# /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/WTDBG/Sebastolobus_alascanus/RagTag/Sebastolobus_altivelis/B-altivelis_xSEB-32/Lifted_Sebastolobus_altivelis.gtf

cat /global/home/users/rohitkolora/RGP/scripts/orthology/illumina_samples_refgen.txt |
 awk -F'\t' '$1!="Sebastes_aleutianus" && $1!="Sebastes_entomelas" && $1!="Sebastes_miniatus" && $1!="Sebastes_pinniger" && $1!="Sebastes_rosaceus" && $1!="Sebastes_umbrosus" && $1!="Sebastolobus_alascanus"' |
 while read line; do
  species=$(echo $line | sed -e 's/ .*//') ;
  sample=$(echo $line | sed -e 's/ /\t/' -e 's/.*\t//' -e 's/ .*//') ;
  refgen=$(echo $line | sed -e 's/.* //') ;
  in_scaff=$(\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/${sample}/ragtag_output/ragtag.scaffolds.fasta | head -n 1) ;

  printf $species" "$in_scaff"\n" ;

done

