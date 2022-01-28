#!/usr/bin/bash

set -e 

work_dir="/global/scratch2/rohitkolora/Rockfish/Transcriptomes/20200521" ;

module load samtools

\ls -d /global/scratch2/rohitkolora/Rockfish/Transcriptomes/20200521/*R1*gz | grep -v -e melano |
    while read fastq1; do
        species_fq=$(basename "$fastq1" | sed -e 's/_R1.fastq.gz$//' -e 's/.*-//' -e 's/rosaceous/rosaceus/' -e 's/miniatius/miniatus/') ;
        outdir=$(basename "$fastq1" | sed -e 's/_R1.fastq.gz$//' -e 's/rosaceous/rosaceus/' -e 's/miniatius/miniatus/') ;
        fastq2=$(echo "$fastq1" | sed -e 's/_R1.fastq.gz$/_R2.fastq.gz/') ;


        \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/referencegenome.fasta | 
            grep -v -e ruber | while read fasta; do 
         
        namen=$(dirname "$fasta" | sed -e 's/\/$//' -e 's/.*\///' )
        species_ref=$(dirname "$fasta" | sed -e 's/\/$//' -e 's/.*\///' -e 's/.*_//') ;
        genome_dir=$(dirname "$fasta" | sed -e 's/$/\/star_index/') ;

        printf " $fasta \t ${species_fq} \t ${species_ref} \n";

        if [[ "${species_fq}" == "${species_ref}" ]]; then
            buty_merged_gtf=$(\ls -d /global/scratch2/rohitkolora/Rockfish/Transcriptomes/gene_coordinates/Final_filt_BRK_*.butyro.merged.gtf | fgrep -e "${species_ref}" | head -n 1) ;
            printf "\t$namen\t$outdir\n" ;
            cd $work_dir/GG_${outdir};

            count_bam=$(find $PWD -type f -name "*.bam" | wc -l) ;
            if [[ "$count_bam" -gt 0 ]] ; then
            bam=$(\ls -d $PWD/*.bam | head -n 1) ;
            count_existing=$(\ls -d $PWD/TPM_butyro_merged/*.out.edit | wc -l) ;
            
            if [[ "$count_existing" -gt 0 ]]; then  # Skip if exists
                continue ;
            fi    

            if [[ -e "$bam" && $bam != "" && ! -f "${bam/.bam/.butyrochr_genes.out}" ]]; then
                printf " Executing TPM calculations\n" ;
                TPMCalculator -g ${buty_merged_gtf} -b ${bam/%.bam/.butyrochr.bam} -p -q 0 -o 25 -e -a ;
                mkdir -p ./TPM_butyro_merged ;
                sed -e "/^FUN_/ s/^/${species_ref}./" ${bam/.bam/.butyrochr_genes.out} >${bam/.bam/.butyrochr_genes.out.edit} ; 
                mv ./*Aligned.sortedByCoord.*butyrochr_genes* ./*Aligned.sortedByCoord.*butyrochr_transcripts* ./TPM_butyro_merged ;
            fi    
            fi    
            cd $work_dir ;
        fi
        done
    done

