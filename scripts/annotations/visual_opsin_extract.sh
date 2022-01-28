#!/usr/bin/bash

set -e

module load samtools bedtools

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/INTERPROSCAN_Fun/interproscan.tsv | 
    while read file1; do 
        species=$(echo $file1 | sed -e 's/.*FREEZE\///' -e 's/\/INTERPRO.*//' -e 's/.*\///'); 
        file2=$(echo $file1 | sed -e "s/INTERPROSCAN_Fun\/interproscan.tsv$/Funannotate\/predict_results\/BRK_${species}.gff3/"); 
        file_prot=$(echo $file1 | sed -e "s/INTERPROSCAN_Fun\/interproscan.tsv$/Funannotate\/predict_results\/Filt_BRK_${species}.proteins.faa/"); 
        file_cds=$(echo $file1 | sed -e "s/INTERPROSCAN_Fun\/interproscan.tsv$/Funannotate\/predict_results\/Filt_BRK_${species}.cds-transcripts.fa/"); 
        file_gff=$(echo $file1 | sed -e "s/INTERPROSCAN_Fun\/interproscan.tsv$/Funannotate\/predict_results\/BRK_${species}.gff3/"); 
        file_ref=$(echo $file1 | sed -e "s/INTERPROSCAN_Fun\/interproscan.tsv$/referencegenome.fasta/");

        cat ~/RGP/scripts/annotations/visualopsin.panther.txt.edit | while read opsindet; do 
            name=$(echo $opsindet | sed -e 's/.* //' -e 's/.*\t//'); 
            opsin=$(echo $opsindet | sed -e 's/ .*//' -e 's/\t.*//' ); 
            #mkdir -p /global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}; 
            #printf ${species}"\t"$name"\n"; 
            #fgrep -w "${opsin}" $file1 | cut -f1 | sort -u | 
            #xargs samtools faidx $file_prot >/global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}/${name}.prot.faa; 
           #fgrep -w "${opsin}" $file1 | cut -f1 | sort -u | 
            #xargs samtools faidx $file_cds >/global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}/${name}.cds.fna;
            #fgrep -w "${opsin}" $file1 | cut -f1 | sort -u |  
            #fgrep -w -f - $file_gff >/global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}/${name}.gff; 
            fgrep -w mRNA /global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}/${name}.gff |
            awk -F'\t' -v name="${name}" 'BEGIN {OFS="\t"} { if($7=="-"){print $1,$4-10001,$5+10001,name"_"NR,$9,$7} else print $1,$4-10001,$5+10001,name"_"NR,$9,$7 }' >/global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}/${name}.bed ;
            bedtools getfasta -fi ${file_ref} -bed /global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}/${name}.bed -s -name >/global/scratch2/rohitkolora/Rockfish/Data/Opsins/annotations/${species}/${name}/${name}.fasta ;
            echo "Done - ${species} : ${name} " ;

        done; 
    done


