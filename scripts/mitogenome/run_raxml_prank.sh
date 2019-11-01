#!/usr/bin/bash

set -e

mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes/raxml
wordir="/global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes"

module load trimal mafft raxml gcc muscle java macse prank
prank -verbose -d=CDS.all.FAS -o=CDS.all.MFA1P.FASTA -DNA -F 
trimal -in CDS.all.MFA1P.FASTA.best.fas -fasta -out CDS.all.MFA2P.FASTA -gt 1 

raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 1 -x 1 -# 100 -s CDS.all.MFA2P.FASTA -n rockfish_mitoP -o Adelosebastes_latens.A_latens_148487

