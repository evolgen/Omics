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
            gtf=$(\ls -d /global/scratch2/rohitkolora/Rockfish/Transcriptomes/gene_coordinates/Final_filt_BRK_*.gtf | fgrep -e "${species_ref}" | fgrep -v '.butyro.gtf' | head -n 1) ;
            printf "\t Extracting Butyro ${gtf}\n";
            cut -f1 /global/scratch2/rohitkolora/Rockfish/Genomes/annotation/rockfish_butyrophilin_genes.vgp_umbrosus_hits.curated.txt | fgrep "${species_ref}" | sed -e 's/.*FUN/FUN/' | sort -u | fgrep -w -f - ${gtf} >${gtf/%.gtf/.butyro.gtf} ;
            printf "\t$namen\t$outdir\n" ;
            cd $work_dir/GG_${outdir};

            count_bam=$(find $PWD -type f -name "*.bam" | wc -l) ;
            if [[ "$count_bam" -gt 0 ]] ; then
            bam=$(\ls -d $PWD/*.bam | head -n 1) ;
            count_existing=$(\ls -d $PWD/TPM_butyro/*.out.edit | wc -l) ;
            
            if [[ "$count_existing" -gt 0 ]]; then  # Skip if exists
                continue ;
            fi    

            if [[ -e "$bam" && $bam != "" && -f "${bam/.bam/_genes.out}" ]]; then
                mkdir -p ./TPM_all ;
                printf " Moving existing TPM data\n" ;
                mv ./*Aligned.sortedByCoord.*out_genes* ./*Aligned.sortedByCoord.*out_transcripts* ./TPM_all ;
            fi    
            if [[ -e "$bam" && $bam != "" && ! -f "${bam/.bam/.butyrochr_genes.out}" ]]; then
                printf " Creating subset bam \n" ;
                cat ${gtf/%.gtf/.butyro.gtf} | cut -f1 | sort -u | fgrep -w -f - ${fasta}.fai | awk -F'\t' '{print $1"\t0\t"$2}' >${gtf/%.gtf/.butyrochr.bed} ;
                samtools view -L ${gtf/%.gtf/.butyrochr.bed} -b -o ${bam/%.bam/.butyrochr.bam} ${bam} ;
                ls -lh ${bam/%.bam/.butyrochr.bam} ;
                samtools index ${bam/%.bam/.butyrochr.bam} ;
                printf " Executing TPM calculations\n" ;
                TPMCalculator -g ${gtf/%.gtf/.butyro.gtf} -b ${bam/%.bam/.butyrochr.bam} -p -q 0 -o 25 -e -a ;
                mkdir -p ./TPM_butyro ;
                sed -e "/^FUN_/ s/^/${species_ref}./" ${bam/.bam/.butyrochr_genes.out} >${bam/.bam/.butyrochr_genes.out.edit} ; 
                mv ./*Aligned.sortedByCoord.*butyrochr_genes* ./*Aligned.sortedByCoord.*butyrochr_transcripts* ./TPM_butyro ;
            fi    
            fi    
            cd $work_dir ;
        fi
        done
    done

