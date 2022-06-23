#!/usr/bin/bash

file1=$@
# ls file1 in new_mod_mafs/Lvir_*.maf | parallel -j 28 run_new_phylop_maf.sh

#rm new_phylop_elements/* new_phylop_scores/* new_phylop_elements_bb/*

file2=$(echo $file1 | sed -e 's/new_mod_mafs/new_phylop_elements/' -e 's/\..*/.element-scores.txt/')
file3=$(echo $file1 | sed -e 's/new_mod_mafs/new_gffs/' -e 's/\..*/.gff/')
file4=$(echo $file1 | sed -e 's/new_mod_mafs/new_phylop_scores/' -e 's/\..*/.scores.wig/')


echo "/scr/bloodymary/rohit/phast/bin/phyloP --features $file3 --method SCORE --mode CONACC Lacerta_init.mod --msa-format MAF $file1 >$file2"
/scr/bloodymary/rohit/phast/bin/phyloP --features $file3 --method SCORE --mode CONACC Lacerta_init.mod --msa-format MAF $file1 >$file2
echo "/scr/bloodymary/rohit/phast/bin/phyloP --wig-scores --method LRT Lacerta_init.mod --msa-format MAF $file1 >$file4"
/scr/bloodymary/rohit/phast/bin/phyloP --wig-scores --method LRT Lacerta_init.mod --msa-format MAF $file1 >$file4



file5=$(echo $file1 | sed -e 's/new_mod_mafs/new_phylop_elements_bb/' -e 's/\..*/.base-scores.wig/')
file6=$(echo $file1 | sed -e 's/new_mod_mafs/new_phylop_elements_bb/' -e 's/\..*/.base-scores.txt/')


echo "/scr/bloodymary/rohit/phast/bin/phyloP --base-by-base --method SCORE --mode CONACC Lacerta_init.mod --msa-format MAF $file1 >$file6"
/scr/bloodymary/rohit/phast/bin/phyloP --base-by-base --method SCORE --mode CONACC Lacerta_init.mod --msa-format MAF $file1 >$file6
echo "/scr/bloodymary/rohit/phast/bin/phyloP --wig-scores --base-by-base --method SCORE --mode CONACC Lacerta_init.mod --msa-format MAF $file1 >$file5"
/scr/bloodymary/rohit/phast/bin/phyloP --wig-scores --base-by-base --method SCORE --mode CONACC Lacerta_init.mod --msa-format MAF $file1 >$file5


echo "$file1	Done"

#; done 


