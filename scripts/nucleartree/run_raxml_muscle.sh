#!/usr/bin/bash

set -e

module load trimal mafft raxml gcc muscle # java macse
muscle -diags -in CDS.FASTA -out CDS.all.MFA1M.FAS
trimal -in CDS.all.MFA1M.FAS -fasta -out CDS.all.MFA2M.FASTA -gt 1


raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 12 -x 12 -# 100 -s CDS.all.MFA2M.FASTA -n rockfish_Muscle -o Adelosebastes_latens

