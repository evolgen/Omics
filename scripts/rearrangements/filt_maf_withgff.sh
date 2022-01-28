#!/usr/bin/bash

set -e -o pipefail

if [ "$#" -ne 3 ]; then
    printf "  Need atleast three arguments - maf, gff, name\n"
    echo "  Run as : bash ~/RGP/scripts/rearrangements/filt_maf_withgff.sh in.maf in.gff name "
    exit 1
fi

maf=$1
gff=$2
name=$3

if [[ ${maf} != *.maf && ${maf} != *.MAF && ${maf} != *.Maf ]]; then
    echo " Please check that maf file name ends with maf"    
    exit 1
fi

if [[ ${gff} != *.gff && ${gff} != *.gff3 && ${gff} != *.GFF && ${gff} != *.GFF3 ]]; then
    echo " Please check that gff file name ends with gff"
    exit 1
fi

if [[ ${name} == "" ]]; then
    echo " Please check that defined name is not empty"
    exit 1
fi

printf "\tOrdering the maf by given species list\n" ;
/global/scratch2/rohitkolora/Software/kentutils/mafOrder ${maf} /global/scratch2/rohitkolora/Rockfish/Data/cactus/species.8ref.list ${name}.order.maf

printf "\tFilter blocks for aleutianus compulsorily\n" ;
/global/scratch2/rohitkolora/Software/kentutils/mafFilter -needComp=Sebastes_aleutianus ${name}.order.maf >${name}.order.keep_aleutianus.maf

printf "\tNeed to remove duplicates blocks in the species - retaining closest to consensus\n" ;
/global/scratch2/rohitkolora/Software/mafTools/bin/mafDuplicateFilter -m ${name}.order.keep_aleutianus.maf >${name}.order.keep_aleutianus.duprem.maf
printf "\tStanding by reference\n" ;
/global/scratch2/rohitkolora/Software/mafTools/bin/mafStrander -m ${name}.order.keep_aleutianus.duprem.maf --seq=Sebastes_aleutianus >${name}.order.keep_aleutianus.duprem.stranded.maf

printf "\tFiltering the maf using gff\n" ;
echo "##maf version=1" >${name}.gff.maf
~/local_modules_sw/phast/1.5/bin/maf_parse -o MAF -g ${gff} ${name}.order.keep_aleutianus.duprem.stranded.maf >>${name}.gff.maf

printf "\tSorting the maf file\n" ;
/global/scratch2/rohitkolora/Software/last-1066/bin/maf-sort ${name}.gff.maf >${name}.final.maf

printf "\n\t OUTPUT is at - ${name}.final.maf \n";


