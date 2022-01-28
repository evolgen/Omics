#!/usr/bin/bash

set -e

file1=$1
workdir=$(dirname "$file1") ;
sample=$(dirname "$file1" | sed -e 's/\/global\/.*\/select_sequence_files\///' -e 's/\/Filter.*//') 

#/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/list_proteins.trimpep.list
#/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/OrthoGroup17684/Filter/aln.pep.out.trimall.fasta

module load samtools trimal gcc ;

printf "  Starting extraction for - ${sample} \n";
            cd ${workdir}/../REVISE ;
            if [[ -f "/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA" ]]; then
                count_revfas=$(grep '^>' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA | wc -l);
                if [[ "$count_revfas" -eq 0 ]]; then
                    rm -f /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA ;
                    exit 0;
                fi    
                trimal -noallgaps -in /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA -out /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA1 ;
                mv /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA1 /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA ;
                count_revfas=$(grep '^>' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA | wc -l);
                if [[ "$count_revfas" -eq 0 ]]; then
                    rm -f /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA ;
                    exit 0;
                fi
                exit 0;
            fi    

            rm -f ../Filter/aln.pep.out.fa.sub ../Filter/aln.cds.orf.aln.fasta.sub ../Filter/aln.pep.out.fa.80.trim ../Filter/aln.cds.orf.aln.fasta.80.trim ;

            if [[ -f "/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/select_filter_trim80/${sample}.FAA" ]]; then
                samtools faidx /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/select_filter_trim80/${sample}.FAA ;

                cat input_species.txt <(cut -f1 /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/select_filter_trim80/${sample}.FAA.fai) |
                    sort | uniq -d |
                    xargs samtools faidx /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/select_filter_trim80/${sample}.FAA >/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA1 ;
                    trimal -noallgaps -in /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA1 -out /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA ;
                rm -f /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/select_filter_trim80/${sample}.FAA.fai /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA1 ;   

                count_newtrim=$(grep '^>' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA | wc -l) ;
                if [[ "$count_newtrim" -lt 10 ]]; then
                    rm -f /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/rerconverge/revise_select_filter_trim80/${sample}.FAA ;
                fi    
            fi

            printf " Trimming sequences\n" ;

            cat input_species.txt | xargs samtools faidx ../Filter/aln.pep.out.fa >../Filter/aln.pep.out.fa.sub ;
            cat input_species.txt | xargs samtools faidx ../Filter/aln.cds.orf.aln.fasta >../Filter/aln.cds.orf.aln.fasta.sub ;

            trimal -fasta -in ../Filter/aln.pep.out.fa.sub -out ../Filter/aln.pep.out.fa.80.trim -resoverlap 0.70 -seqoverlap 80 ;
            trimal -fasta -in ../Filter/aln.cds.orf.aln.fasta.sub -out ../Filter/aln.cds.orf.aln.fasta.80.trim -resoverlap 0.70 -seqoverlap 80 ;

            printf " Indexing trimmed sequences\n" ;

            if [[ ! -f "../Filter/aln.pep.out.fa.80.trim" ]]; then
                echo "$sample" >group_trimming_prob.txt ;
                exit 1;
            fi

            if [[ ! -f "../Filter/aln.cds.orf.aln.fasta.80.trim" ]]; then
                echo "$sample" >group_trimming_prob.txt ;
                exit 1;
            fi    

            samtools faidx ../Filter/aln.pep.out.fa.80.trim && samtools faidx ../Filter/aln.cds.orf.aln.fasta.80.trim ;
            cat input_species.txt <(cut -f1 ../Filter/aln.pep.out.fa.80.trim.fai) <(cut -f1 ../Filter/aln.cds.orf.aln.fasta.80.trim.fai) |
                sort | uniq -d | grep -w 3 | sed -e 's/.*Seb/Seb/' >select_species.txt ;

            count_species=$(cat select_species.txt | wc -l) ;
            printf "    : extraction\n" ;
            if [[ "$count_species" -gt 10 ]]; then
                cat select_species.txt | xargs samtools faidx ../Filter/aln.pep.out.fa.80.trim >80tr.PEP.fas ;
                cat select_species.txt | xargs samtools faidx ../Filter/aln.cds.orf.aln.fasta.80.trim >80tr.CDS.fas ;
            else
                exit 1;
            fi

            samtools faidx 80tr.PEP.fas && samtools faidx 80tr.CDS.fas ;
            pep_size=$(cut -f2 80tr.PEP.fas.fai | sort -n | head -n 1) ;
            cds_size=$(cut -f2 80tr.CDS.fas.fai | sort -n | head -n 1) ;

            printf "    : validation\n" ;
            if [[ "$pep_size" -lt 40 && "$cds_size" -lt 120 ]]; then
                rm -f 80tr.PEP.fas 80tr.PEP.fas.fai 80tr.CDS.fas 80tr.CDS.fas.fai ;
            fi

