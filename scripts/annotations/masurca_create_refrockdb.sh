#!/usr/bin/bash

set -e -o pipefail

module load samtools bedtools

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/assembly/masurca/Funannotate/predict_results/Filt_BRK_*proteins.faa | #head -n 1 |
    grep -v S-ruberrimus_SEB-74 | while read file; do 
        echo ${file} ; printf "\n\n" ;

        species=$(basename "${file}" | sed -e 's/.*BRK_//' -e 's/\..*//') ;
        outfile=$(echo $file | sed -e 's/\/Filt_/\/Final_filt_/') ;
        cdsfile=$(echo $file | sed -e 's/\/Filt_/\//' -e 's/\.proteins.faa$/.cds-transcripts.fa/') ;
        outcdsfile=$(echo $file | sed -e 's/\/Filt_/\/Final_filt_/' -e 's/\.proteins.faa$/.cds-transcripts.fa/') ;
        gff3=$(echo $file | sed -e 's/\/Filt_/\//' -e 's/.proteins.faa/.gff3/')
        outgff=$(echo $file | sed -e 's/\/Filt_/\/Final_filt_/' -e 's/.proteins.faa/.gff3/') ;

        diamond blastp -d /global/scratch2/rohitkolora/databases/rockfish/refgen_proteins_fun -q ${file} --threads $(eval nproc) --outfmt 6 --top 0.1 --unal 1 --index-chunks 1 -o $(dirname "${file}")/simhits_otherrockref.txt ;
        printf "\tBlast done\n";

        awk '$2!="*" && !x[$1]++' $(dirname "${file}")/simhits_otherrockref.txt | 
            cut -f1 | 
            xargs samtools faidx ${file} >${outfile}.out.pep ;
        samtools faidx ${outfile}.out.pep ;    
        printf "\tPep extract done\n";

        cut -f1 ${outfile}.out.pep.fai >${outfile}.out.pepnames ;
        grep -w -e 'mRNA' ${gff3} |
            grep -v _assoc | awk -F'\t' '{print $NF}' |
            sed -e 's/ID=//' -e 's/;.*//' | sort -u |
            fgrep -w -f - ${outfile}.out.pepnames >${outfile}.out.names ;
        printf "\tName extract based on mRNA done\n";    

        cat ${outfile}.out.names |
            xargs samtools faidx ${outfile}.out.pep >${outfile} ;
        printf "\tSequence extract done\n" ;    
        samtools faidx ${outfile} ;

        cut -f1 ${outfile}.fai |
            xargs samtools faidx ${cdsfile} >${outcdsfile} ;
        printf "\tRunning gff extract\n"; 
        cut -f1 ${outfile}.fai | 
            awk '{print $1; gsub("-T[0-9]*", "", $1); print $1}' |
            fgrep -w -f - ${gff3} | sort -u |
            sort -k1,1V -k4,4n -k5,5nr -k3,3V >${outgff} ;
        printf "\t$species\t" ; wc -l ${outfile}.fai ;
        echo ;

    done

