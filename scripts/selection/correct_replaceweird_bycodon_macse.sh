#!/usr/bin/bash

set -e

module load java macse trimal gcc

probfile=$@

#count=0; find /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/nucleargenome/SCOs_76/sequence_files/ -type f -name "trim.NT.seqerr" | while read fasta; do prob=$(dirname $fasta | sed -e 's/$/\/trim.NT.seqerr/') if [[ -e "" ]]; then then count=$((count+1)); trimal -gt 1 -in $fasta -out -nexus ; echo $fasta"    "$count; fi; done

fasta=$(dirname $probfile | sed -e 's/$/\/seq.fasta/')

echo $fasta
java -Xmx30000m -jar /global/home/users/rohitkolora/local_modules_sw/macse/macse_v2.03.jar -prog alignSequences -seq $fasta -local_realign_init 1 -local_realign_dec 1 -seq $fasta -out_NT ${fasta/seq.fasta/aln.NT.fasta} -out_AA ${fasta/seq.fasta/aln.AA.fasta}

trimal -gt 1 -in ${fasta/seq.fasta/aln.NT.fasta} -out ${fasta/seq.fasta/trim.NT.fasta} -fasta

perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${fasta/seq.fasta/trim.NT.fasta} | sed -e '/^>/! s/NNN$//' >${fasta/seq.fasta/trim.NT.fasta.edit} 
trimal -gt 1 -in ${fasta/seq.fasta/trim.NT.fasta.edit} -out ${fasta/seq.fasta/trim.NT.nex} -nexus 1>${fasta/seq.fasta/trim.NT.log}   

value=$(grep -m 1 '#NEXUS' ${fasta/seq.fasta/trim.NT.nex} | wc -l)
if [[ "${value}" -gt 0 ]]; then
    rm -f $probfile
fi    

