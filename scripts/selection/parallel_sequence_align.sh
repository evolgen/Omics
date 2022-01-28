#!/usr/bin/bash

set -e -o pipefail

module load samtools trimal mafft raxml gcc muscle java macse;
#conda activate hyphy;

mkdir -p $PWD/sequence_files #$PWD/hyphy_busted;
#workingdir="/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome"

infile=$@
genefile=$(basename "${infile}" | sed -e 's/.fasta$//' -e 's/.FASTA$//' -e 's/.fas$//' -e 's/.fas$//' -e 's/.faa$//' -e 's/.FAA$//' -e 's/.fna$//' -e 's//.FNA$/' -e 's/.fa$//' -e 's/.FA$//')

printf "" >$PWD/sequence_files/${genefile}.fasta ;
echo "${genefile}" ;

grep '^>' $infile | grep -e 'Sebastes_' |
    sed -e 's/^>//' -e 's/ .*//' |
    sort -u |
    while read speciesname; do
        printf ">${speciesname}\n" >> $PWD/sequence_files/${genefile}.fasta ;
        samtools faidx $infile ${speciesname} |
            grep -v '^>' | 
            sed -e '/^$/d' -e '$s/TAA$/NNN/' -e '$s/TAG$/NNN/' -e '$s/TGA$/NNN/' -e '$s/taa$/NNN/' -e '$s/tag$/NNN/' -e '$s/tga$/NNN/' |
            tr "\n" " " | sed -e 's/ //g' -e 's/$/\n/' |
            fold -w 3 |
            awk 'length($0) == 3' |
            tr "\n" " " | sed -e 's/ //g' -e 's/$/\n/' >>$PWD/sequence_files/${genefile}.fasta;
    done          

echo;

muscle -diags -in $PWD/sequence_files/${genefile}.fasta -out $PWD/sequence_files/${genefile}.aln;
trimal -gt 1 -in $PWD/sequence_files/${genefile}.aln -nexus -out $PWD/sequence_files/${genefile}.nex; 
        

