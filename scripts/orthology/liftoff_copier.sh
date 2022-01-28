#!/usr/bin/sh

set -e -o pipefail

#ls /global/home/users/rohitkolora/RGP/scripts/orthology/illumina_samples_refgen.txt
# /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/WTDBG/Sebastolobus_alascanus/RagTag/Sebastolobus_altivelis/B-altivelis_xSEB-32/Lifted_Sebastolobus_altivelis.gtf

mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/input_data ;
cat /global/home/users/rohitkolora/RGP/scripts/orthology/illumina_samples_refgen.txt |
 awk -F'\t' '$1!="Sebastes_aleutianus" && $1!="Sebastes_entomelas" && $1!="Sebastes_miniatus" && $1!="Sebastes_pinniger" && $1!="Sebastes_rosaceus" && $1!="Sebastes_umbrosus" && $1!="Sebastolobus_alascanus"' |
 while read line; do
  species=$(echo $line | sed -e 's/ .*//') ;
  sample=$(echo $line | sed -e 's/ /\t/' -e 's/.*\t//' -e 's/ .*//') ;
  refgen=$(echo $line | sed -e 's/.* //') ;
  in_gff=$(\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/${sample}/${species}.gff | head -n 1) ;
  in_prot=$(\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/${sample}/${species}.faa | head -n 1) ;
  in_cds=$(\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/${sample}/${species}.fna | head -n 1) ;
#    in_gff="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/${sample}/${species}.gff" ;
#    in_prot="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/${sample}/${species}.faa" ;
#    in_cds="/global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/${sample}/${species}.fna" ;
#if [[ ! -e "$in_gff" ]]; then  
  sed -e "/^#/! s/^/${species}./" $in_gff |
    awk -v species="$species" -F'\t' 'BEGIN{OFS="\t"} {gsub("ID=","ID="species".",$9); print $0}' >/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/input_data/${species}.gff ;
  sed -e "/^>/ s/^>/>${species}./" $in_prot >/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/input_data/${species}.faa ;
  sed -e "/^>/ s/^>/>${species}./" $in_cds >/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/input_data/${species}.fna ;
  printf "  Done with\t$species\t\t$sample\t\tBase : $refgen\n" ;
#fi
 done    

#\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/proteinortho_finale/input_data/*faa /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/proteinortho_finale/input_data/*gff | grep -e Sebastes_aleutianus -e Sebastes_entomelas -e Sebastes_miniatus -e Sebastes_pinniger -e Sebastes_rosaceus -e Sebastes_umbrosus -e Sebastolobus_alascanus | while read file; do cp $file /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/input_data; done

