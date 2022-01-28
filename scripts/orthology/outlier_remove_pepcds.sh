#!/usr/bin/bash

set -e

#cat /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/list_proteins.trimpep.list | parallel -j 1 bash ~/RGP/scripts/orthology/outlier_from_pepcds.sh
#/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/OrthoGroup19105/Filter/aln.pep.out.trimall.fasta

in_pep_fa=$1

in_cds_fa=$(echo $in_pep_fa | sed -e 's/aln.pep.out.trimall.fasta/aln.cds.orf.aln.fasta/') ;

orthogroup=$(dirname "$in_pep_fa" | sed -e 's/\/Filter$//' -e 's/.*\/select_sequence_files\///') ;

out_dir=$(dirname "$in_pep_fa" | sed -e 's/\/Filter$/\/REVISE/') ;

cd $out_dir ;

echo $orthogroup ;
module load samtools ;

samtools faidx ../Filter/aln.cds.out.80.trim.fas && samtools faidx ../Filter/aln.pep.out.80.trim.fas ;

cat input_species.txt | xargs samtools faidx ../Filter/aln.cds.out.80.trim.fas >${out_dir}/ALN.CDS.80trim.fasta ;
cat input_species.txt | xargs samtools faidx ../Filter/aln.pep.out.80.trim.fas >${out_dir}/ALN.PEP.80trim.fasta ;

samtools faidx ${out_dir}/ALN.CDS.80trim.fasta && samtools faidx ${out_dir}/ALN.PEP.80trim.fasta ;   

