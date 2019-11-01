#!/usr/bin/bash

set -e

mkdir -p /global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes/raxml
wordir="/global/scratch2/rohitkolora/Rockfish/Genomes/mitogenomes"

#touch ${wordir}/CDS.all.FAS
##\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/MITOS2/CDS/ALL.FAS |
##    xargs cat | 
##    sed -e 's/>/\n>/' | 
##    sed -e '/^$/d' >${wordir}/CDS.all.FAS

module load trimal mafft raxml gcc muscle # java macse
mafft --thread 32 --maxiterate 1000 --localpair ${wordir}/CDS.all.FAS >${wordir}/CDS.all.MFA
#java -Xmx300g -jar /global/home/users/rohitkolora/local_modules_sw/macse/macse.jar -i CDS.all.FAS -d 2 -o rockfish_mito -o_dna CDS.all.MFA
#prank -d=CDS.all.FAS -o=CDS.all.MFA1P.FASTA -DNA -F && trimal -in CDS.all.MFA1P.FASTA.best.fas -fasta -out CDS.all.MFA2P.FASTA -gt 1
trimal -in ${wordir}/CDS.all.MFA -fasta -out ${wordir}/CDS.all.MFA1.FAS -gt 1
muscle -diags -in ${wordir}/CDS.all.MFA1.FAS -out ${wordir}/CDS.all.MFA2.FAS
trimal -in ${wordir}/CDS.all.MFA2.FAS -fasta -out ${wordir}/CDS.all.MFA2D.FASTA -gt 1


grep '^>' ${wordir}/CDS.all.MFA.FASTA | 
    sed -e 's/^>//' | 
    tr "\n" "," | sed -e 's/,$/\n/'

number=$(date '+%d/%m/%Y/%H:%M:%S');

raxmlHPC-SSE3 -T 32 -f a -m GTRGAMMA -p 1234 -x 1234 -# 100 -s ${wordir}/CDS.all.MFA2D.FASTA -n rockfish_mitoD -o Adelosebastes_latens.A_latens_148487

