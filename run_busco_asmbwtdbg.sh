#!/usr/bin/bash

set -e

filepattern=$1
outname=$2

module load busco/3.1 augustus/2.5.5 hmmer blast/2.2.26 gcc
AUGUSTUS_CONFIG_PATH="/global/home/users/rohitkolora/local_modules_sw/busco/3.1//../../augustus/2.5.5/config/"

if [ "$#" -ne 2 ]; then
    printf "Illegal number of parameters:\t     $#      \n";
    printf "\tPlease RUN as : sh run_busco_asmbwtdbg.sh FilePattern OutPattern\n\tNOTE: Only fasta files with *.fa extension allowed";
    printf "\t\tExample:        sh run_busco_asmbwtdbg.sh /PATH/TO/FILES/PATTERN OUTFILEPATTERN\n";                         #sh extract_pacbiorun_tar.sh $PWD/wtdbg2_007.V run_busco 
    exit;
fi

FASTA_PATHING=$1"*.fa"
echo "Running with files as $FASTA_PATHING"

OUTPATTERN=$2

checking=$(ls $FASTA_PATHING | grep -e '\.fa$' -e '\.fasta$' | wc -l)
if [ "$checking" -eq 0 ]; then
    printf "FASTA files with pattern $FASTA_PATHING dont exist\n\tPLEASE CHECK for $1*.fa\n";
    exit
fi


for file in $FASTA_PATHING; do
    out=$(echo $file | sed -e 's/.*\///' -e 's/\.V/_V/' -e 's/\.fa.*/_acti/' -e 's/ctg//' -e 's/\.//g'); 
    python /global/home/users/rohitkolora/local_modules_sw/busco/3.1//scripts/run_BUSCO.py -i $file -c 32 -o ${OUTPATTERN}_${out} -m geno -sp zebrafish -l /global/scratch/rohitkolora/databases/busco/actinopterygii_odb9 ; 
done


