#!/usr/bin/bash

set -e

mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes/raxml
wordir="/global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes"

module load trimal mafft raxml gcc muscle # java macse
muscle -diags -in ${wordir}/CDS.all.FAS -out ${wordir}/CDS.all.MFA1M.FAS
trimal -in ${wordir}/CDS.all.MFA1M.FAS -fasta -out ${wordir}/CDS.all.MFA2M.FASTA -gt 1

number=$(date '+%d/%m/%Y/%H:%M:%S');

raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 12 -x 12 -# 100 -s ${wordir}/CDS.all.MFA2M.FASTA -n rockfish_mitoM -o Adelosebastes_latens.A_latens_148487

