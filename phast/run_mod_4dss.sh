#!/bin/bash

set -e

file1=$@

# for file1 in ./mafs/Lvir_*.maf; do 
infile=$(echo $file1 | sed -e 's/mod_mafs/seqs/' -e 's/\.maf/.fa/')
file2=$(echo $file1 | sed -e 's/mod_mafs/check_gffs/' -e 's/\.maf/.gff/')
file3=$(echo $file2 | sed -e 's/check_gffs/check_codons_ss/' -e 's/\.gff/.codons.ss/')
file4=$(echo $file3 | sed -e 's/check_codons_ss/check_sites_ss/' -e 's/\.ss/.sites.ss/')


if [ -f $file2 ]
then

###msa_view $infile --4d --features $file2 > $file3

# msa_view $file1 --in-format MAF --out-format SS --4d --features $file2 > $file3		# For codon
# msa_view $file3 --in-format SS --out-format SS --tuple-size 1 > $file4			# for sites


fi

#else
#	echo "$file2 doesn't exist."
#	echo $file1 >>missing_list.txt
#fi



#; done


