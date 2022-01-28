#!/usr/bin/sh

set -e -o pipefail

if [ "$#" -ne 3 ]; then
          echo "Input arguments are syntfile chrnameinfo outfile";
          exit 1;
fi

syntfile1=$1
chrnameinfo=$2
outfile=$3

syntfile_subset="${syntfile1/%.out/.syntenies}"
awk '$11=="SYN"' $syntfile1 >$syntfile_subset ;
printf "  Created subset syntenic file\n";

refreplacerfile="${chrnameinfo}.refnames"
sed '/###REF-names###/,/###QRY-names###/!d' ${chrnameinfo} |
    awk '$1!~"^#" && NF==2' | 
    sed -e 's/>//g' >${refreplacerfile} ;
printf "  Created REF replacer file\n";

qryreplacerfile="${chrnameinfo}.qrynames"
sed -n '/^###QRY-names###$/,$p' ${chrnameinfo} |
    awk '$1!~"^#" && NF==2' |
    sed -e 's/>//g' >${qryreplacerfile} ;
printf "  Created QRY replacer file\n";

count_syntfile=$(cat $syntfile_subset | awk '{print NF}' | sort -u);
count_refreplacerfile=$(cat $refreplacerfile | awk '{print NF}' | sort -u);
count_qryreplacerfile=$(cat $qryreplacerfile | awk '{print NF}' | sort -u);

printf "  Counted lines for inputs\n" ;

if [ "$count_syntfile" -ne 12 ] || [ "$count_refreplacerfile" -ne 2 ] || [ "$count_qryreplacerfile" -ne 2 ] ; then
    echo "Input file errors for syntfile, refreplacerfile OR qryreplacerfile";
    exit 1;
fi

#source ~/conda.sh && conda activate R
printf "  Executing R script for replacement\n"
Rscript ~/RGP/scripts/rearrangements/replace_synt_withlist.R $syntfile_subset $refreplacerfile $qryreplacerfile $outfile

