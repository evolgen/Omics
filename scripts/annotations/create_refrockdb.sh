#!/usr/bin/bash

set -e -o pipefail

module load samtools bedtools


\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/Funannotate/predict_results/filt_BRK_*proteins.fa |
    grep -v ruberrimus | while read file; do 
        echo ${file} ; 
        folder=$(dirname "${file}" | sed -e 's/\/predict_results.*//') ; 
        species=$(basename "${file}" | sed -e 's/.*BRK_//' -e 's/\..*//') ;
        outfile=$(echo $file | sed -e 's/\/filt_/\/Final_filt_/') ;
        cdsfile=$(echo $file | sed -e 's/\/filt_/\//' -e 's/\.proteins.fa$/.cds-transcripts.fa/') ;
        outcdsfile=$(echo $file | sed -e 's/\/filt_/\/Final_filt_/' -e 's/\.proteins.fa$/.cds-transcripts.fa/') ;
        gff3=$(echo $file | sed -e 's/\/filt_/\//' -e 's/.proteins.fa/.gff3/')
        outgff=$(echo $file | sed -e 's/\/filt_/\/Final_filt_/' -e 's/.proteins.fa/.gff3/') ;
        \ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/Funannotate/predict_results/filt_BRK_*proteins.fa | 
            grep -v ruberrimus | grep -v "${file}" | 
            xargs cat >${folder}/other_ref-rockfish.prot.faa ;
        cd ${folder}/ ;   
        diamond makedb --in other_ref-rockfish.prot.faa -d other_ref-rockfish_prot_fun ;  
        diamond blastp -d other_ref-rockfish_prot_fun -q ${file} --threads $(eval nproc) --outfmt 6 --top 0.1 --unal 1 --index-chunks 1 -o $(dirname "${file}")/simhits_otherrockref.txt ;
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
        cp ${outfile} ${outfile}a ;
        cut -f1 ${outfile}.fai |
            xargs samtools faidx ${cdsfile} >${outcdsfile} ;
        samtools faidx ${outfile}a ;    

        samtools faidx ${outfile} ;
        printf "\tRunning gff extract\n"; 
        cut -f1 ${outfile}.fai | 
            awk '{print $1; gsub("-T[0-9]*", "", $1); print $1}' |
            fgrep -w -f - ${gff3} | sort -u |
            sort -k1,1V -k4,4n -k5,5nr -k3,3V >${outgff} ;
    done

printf "\tNow executing all reference genomes\n" ;

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/Funannotate/predict_results/Final_filt_BRK_*proteins.fa |
    grep -v -e ruberrimus | xargs cat >/global/scratch2/rohitkolora/databases/rockfish/refgen_proteins.faa ;
diamond makedb --in /global/scratch2/rohitkolora/databases/rockfish/refgen_proteins.faa -d /global/scratch2/rohitkolora/databases/rockfish/refgen_proteins_fun ;



