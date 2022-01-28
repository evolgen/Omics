#!/usr/bin/bash

#set -e -o pipefail

module load samtools ;

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/RagTag/*/*/Lifted_*.gtf |
    while read gtf; do
        echo "${gtf}" ;
        prot_in=$(echo ${gtf} | sed -e 's/.gtf$/.proteins.faa/') ;
        cds_in=$(echo ${gtf} | sed -e 's/.gtf$/.cds-transcripts.fa/') ;
        prot_seq=$(echo ${prot_in} | sed -e 's/\/Lifted_/\//' -e 's/.proteins.faa$/.faa/') ;
        gff=$(echo ${gtf} | sed -e 's/\/Lifted_/\//' -e 's/.gtf$/.gff/') ;
        cds_seq=$(echo ${cds_in} | sed -e 's/\/Lifted_/\//' -e 's/.cds-transcripts.fa$/.fna/') ;
        seq_names=$(echo ${gtf} | sed -e 's/\/Lifted_/\//' -e 's/.gtf$/.names/') ;
#    if [[ ! -e "${gff}" ]]; then
            printf "\tExtracting names\t";
        #awk -F'\t' '$3=="transcript" {print $NF}' ${gtf} |
        #    sed -e 's/transcript_id "//' -e 's/";.*//' |
        #    sort -u >${seq_names} ;
        samtools faidx ${prot_in} ;  
        awk -F'\t' '$2>=30 {print $1}' ${prot_in}.fai | 
            sort -u >${seq_names} ;
        printf "\tExtracting names\t";
        cat ${seq_names} | 
            xargs samtools faidx ${prot_in} >${prot_seq} ;
        samtools faidx ${prot_seq} ;
        printf "Extracting cds\t";
        samtools faidx ${cds_in} ;
        cat ${seq_names} | 
            xargs samtools faidx ${cds_in} >${cds_seq} ;
        samtools faidx ${cds_seq} ;
        printf "Extracting gff\t";
        cat ${seq_names} |
            fgrep -w -f - ${gtf} |
            awk -F'\t' '$3=="CDS"' |
            sed -e 's/transcript_id "/ID=/' -e 's/";.*/;/' >${gff} ;
        #awk -F'\t' '$3=="CDS"' ${gtf} |
        #    sed -e 's/transcript_id "/ID=/' -e 's/";.*/;/' |
        #    fgrep -w -f <(cut -f1 ${prot_seq}.fai | sort -u) >${gff} ;
        printf "\n\tProcessed : ${gtf}\n\n";    
#    fi
#        printf "  $gtf\t$prot_in\t$prot_seq\t$cds_in\t$cds_seq\t$gff\n";
    done 


