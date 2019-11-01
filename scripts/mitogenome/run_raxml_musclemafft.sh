#!/usr/bin/bash

set -e

mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes/raxml
wordir="/global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes"

#touch ${wordir}/CDS.all.FAS

module load trimal mafft raxml gcc muscle # java macse
muscle -diags -in ${wordir}/CDS.all.FAS -out ${wordir}/CDS.all.MFA1R.FAS
trimal -in ${wordir}/CDS.all.MFA1R.FAS -fasta -out ${wordir}/CDS.all.MFA1.FASTA -gt 1
mafft --thread 32 --maxiterate 1000 --localpair ${wordir}/CDS.all.MFA1.FASTA >${wordir}/CDS.all.2R.MFA
trimal -in ${wordir}/CDS.all.2R.MFA -fasta -out ${wordir}/CDS.all.MFA2R.FASTA -gt 1

raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 123 -x 123 -# 100 -s ${wordir}/CDS.all.MFA2R.FASTA -n rockfish_mitoR -o Adelosebastes_latens.A_latens_148487

