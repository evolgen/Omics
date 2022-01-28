#!/usr/bin/bash

# /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/RagTag/Sebastes_aurora/S-aurora_KUI448/Lifted_Sebastes_aurora.proteins.faa
# \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/FALCON/Sebastes_aleutianus/RagTag/*/*/Lifted_*.proteins.faa | wc -l
# /global/home/users/rohitkolora/RGP/scripts/orthology/illumina_samples_refgen.txt

module load samtools ;
mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted ;

awk -F'\t' '$1==$3 {print $1"\t"$3}' /global/home/users/rohitkolora/RGP/scripts/orthology/illumina_samples_refgen.txt |
  while read line1; do
    species=$(echo $line1 | sed -e 's/ .*//') ;
    refgen=$(echo $line1 | sed -e 's/.* //') ;

    mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/ ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${species}/Funannotate/predict_results/Final_filt_BRK_${species}.proteins.faa | sed -e "/^>/ s/>FUN_/>${species}.FUN_/" >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.faa ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${species}/Funannotate/predict_results/Final_filt_BRK_${species}.cds-transcripts.fa | sed -e "/^>/ s/>FUN_/>${species}.FUN_/" >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.fna ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${species}/Funannotate/predict_results/Final_filt_BRK_${species}.nameedit.gff3 | sed -e "s/FUN_/${species}.FUN_/g" >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.gff3 ;

      \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${species}/Pairwise/Liftoff/*_*_*/Lifted_*.proteins.faa |
      grep -v "${refgen}" |
      while read in_prots; do
          in_cds=$(echo $in_prots | sed -e 's/.proteins.faa$/.cds-transcripts.fa/') ;
          in_gff3=$(echo $in_prots | sed -e 's/.proteins.faa$/.gff3/') ;
          in_gtf=$(echo $in_prots | sed -e 's/.proteins.faa$/.gtf/') ;
          refgen2=$(basename "$in_prots" | sed -e 's/.proteins.faa$//' -e 's/Lifted_//') ;

          cat $in_prots | sed -e "/^>/ s/>FUN_/>${species}.${refgen2}.FUN_/" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.faa ;
          cat $in_cds | sed -e "/^>/ s/>FUN_/>${species}.${refgen2}.FUN_/" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.fna ;
          cat $in_gff3 | sed -e "s/FUN_/${species}.${refgen2}.FUN_/g" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.gff3 ;
          cat $in_gtf | sed -e "s/FUN_/${species}.${refgen2}.FUN_/g" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.gtf ;
      done
          samtools faidx /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.faa ;
          samtools faidx /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.fna ;
          echo $species ;
  done    

awk -F'\t' '$1!=$3 {print $1"\t"$3}' /global/home/users/rohitkolora/RGP/scripts/orthology/illumina_samples_refgen.txt |
  while read line1; do
    species=$(echo $line1 | sed -e 's/ .*//') ;
    refgen=$(echo $line1 | sed -e 's/.* //') ;

    mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/ ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/*/Lifted_${species}.proteins.faa | sed -e "/^>/ s/>FUN_/>${species}.FUN_/" >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.faa ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/*/Lifted_${species}.cds-transcripts.fa | sed -e "/^>/ s/>FUN_/>${species}.FUN_/" >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.fna ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/*/Lifted_${species}.gff3 | sed -e "s/FUN_/${species}.FUN_/g" >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.gff3 ;
    cat /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/${refgen}/RagTag/${species}/*/Lifted_${species}.gtf | sed -e "s/FUN_/${species}.FUN_/g" >/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.gtf ;

    \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/RagTag/${species}/*/Lifted_${species}.proteins.faa |
      grep -v "${refgen}" |
      while read in_prots; do
          in_cds=$(echo $in_prots | sed -e 's/.proteins.faa$/.cds-transcripts.fa/') ;
          in_gff3=$(echo $in_prots | sed -e 's/.proteins.faa$/.gff3/') ;
          in_gtf=$(echo $in_prots | sed -e 's/.proteins.faa$/.gtf/') ;
          refgen2=$(dirname "$in_prots" | sed -e 's/.*\/FALCON\///' -e 's/.*\/WTDBG\///' -e 's/\/.*//') ;

          cat $in_prots | sed -e "/^>/ s/>FUN_/>${species}.${refgen2}.FUN_/" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.faa ;
          cat $in_cds | sed -e "/^>/ s/>FUN_/>${species}.${refgen2}.FUN_/" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.fna ;
          cat $in_gff3 | sed -e "s/FUN_/${species}.${refgen2}.FUN_/g" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.gff3 ;
          cat $in_gtf | sed -e "s/FUN_/${species}.${refgen2}.FUN_/g" >>/global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.gtf ;
      done
          samtools faidx /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.faa ;
          samtools faidx /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/All_lifted/${species}/${species}.fna ;
          echo $species ;
  done    
  
