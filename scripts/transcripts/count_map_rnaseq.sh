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

        if [[ "${species_fq}" == "${species_ref}" ]]; then
            gtf=$(\ls -d /global/scratch2/rohitkolora/Rockfish/Transcriptomes/gene_coordinates/Final_filt_BRK_*.gtf | fgrep -e "${species_ref}" | head -n 1) ;
            printf $namen"\t"$outdir"\n" ;
            cd $work_dir/GG_${outdir};

            count_bam=$(find $PWD -type f -name "*.bam" | wc -l) ;
            if [[ "$count_bam" -gt 0 ]] ; then
            bam=$(\ls -d $PWD/*.bam | head -n 1) ;
            if [[ -e "$bam" && $bam != "" && ! -f "${bam/.bam/_genes.out}" ]]; then
                samtools index $bam ;
                TPMCalculator -g ${gtf} -b ${bam} -p -q 0 -o 25 -e -a ;
            fi    
            else
                STAR --runMode alignReads --runThreadN $(eval nproc) --genomeDir $genome_dir --readFilesIn $fastq1 $fastq2 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMprimaryFlag AllBestScore --readFilesCommand zcat --twopassMode Basic --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 # --outFileNamePrefix star_${namen}__${outdir} ;
                samtools index Aligned.sortedByCoord.out.bam ;
                bam=$(\ls -d $PWD/*.bam | head -n 1) ;
                TPMCalculator -g ${gtf} -b ${bam} -p -q 0 -o 25 -e -a ;
            fi    
            cd $work_dir ;
        fi
        done
    done

