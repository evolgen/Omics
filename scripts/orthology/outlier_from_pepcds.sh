#!/usr/bin/bash

set -e

#cat /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/list_proteins.trimpep.list | parallel -j 1 bash ~/RGP/scripts/orthology/outlier_from_pepcds.sh
#/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/OrthoGroup19105/Filter/aln.pep.out.trimall.fasta

in_pep_fa=$1

in_cds_fa=$(echo $in_pep_fa | sed -e 's/aln.pep.out.trimall.fasta/aln.cds.orf.aln.fasta/') ;

orthogroup=$(dirname "$in_pep_fa" | sed -e 's/\/Filter$//' -e 's/.*\/select_sequence_files\///') ;

out_dir=$(dirname "$in_pep_fa" | sed -e 's/\/Filter$/\/REVISE/') ;

mkdir -p $out_dir && cd $out_dir ;
echo $orthogroup ;
module load samtools ;

/global/scratch2/rohitkolora/Software/OD-Seq/OD-seq -i ${in_pep_fa} -t 1 --boot-rep 1000 -r ${out_dir}/outlier_pep.txt -c ${out_dir}/Aln.PEP.trimall.FAS ;

/global/scratch2/rohitkolora/Software/OD-Seq/OD-seq -i ${in_cds_fa} -t 1 --boot-rep 1000 -r ${out_dir}/outlier_cds.txt -c ${out_dir}/Aln.CDS.trimall.FAS ; 

samtools faidx ${out_dir}/Aln.PEP.trimall.FAS && samtools faidx ${out_dir}/Aln.CDS.trimall.FAS ;
count_cds=$(cat ${out_dir}/Aln.PEP.trimall.FAS | wc -l) ;
count_pep=$(cat ${out_dir}/Aln.CDS.trimall.FAS | wc -l) ;

if [[ "${count_cds}" == "${count_pep}" ]]; then
    echo "$orthogroup" >good_group.txt ;
else
    echo "$orthogroup" >bad_group.txt ;
fi    

cat ${out_dir}/Aln.PEP.trimall.FAS.fai ${out_dir}/Aln.CDS.trimall.FAS.fai | 
    cut -f1 | sort | uniq -d >input_species.txt ;

cat input_species.txt | xargs samtools faidx ${out_dir}/Aln.CDS.trimall.FAS >${out_dir}/ALN.CDS.trimall.FAS ;
cat input_species.txt | xargs samtools faidx ${out_dir}/Aln.PEP.trimall.FAS >${out_dir}/ALN.PEP.trimall.FAS ;

samtools faidx ${out_dir}/ALN.CDS.trimall.FAS && samtools faidx ${out_dir}/ALN.PEP.trimall.FAS ;   

