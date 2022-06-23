#!/usr/bin/sh

set -e 

 /scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/mafSplit -byTarget -useFullSequenceName lacvir1.bed ./new_mafs/ ./Lacerta_multiple_lacvir1_roast.MAF
 printf "Broke MAF files per chromsome\n"

 for file1 in ./new_mafs/Lvir_*.maf; do file2=$(echo $file1 | sed 's/new_mafs/new_mod_mafs/'); /scr/bloodymary/rohit/phast/bin/maf_parse -o MAF -O lacvir1,lacbil1,anoCar2 $file1 >$file2; done
printf "Parsed MAF files per chromsome\n"


# grep -v '^#' viridis_CDS_new.gff | awk '{print >> "./new_gffs/"$1".gff"}'; for file in new_gffs/Lvir_*; do sed -i 's/^/lacvir1./' $file; done
# mkdir -p other_gffs
# grep -v '^#' new_CDS_viridis.gff | awk '{print >> "./other_gffs/"$1".gff"}'; for file in other_gffs/Lvir_*; do sed -i 's/^/lacvir1./' $file; done
printf "Created GFF files per chromsome\n"


mkdir -p new_codons_ss new_sites_ss

 ls ./new_mod_mafs/Lvir_*.maf | parallel -j 30 sh run_new_con_4dss.sh
printf "4dss files created\n"

/scr/bloodymary/rohit/phast/bin/maf_parse -o SS -O lacvir1,lacbil1,anoCar2 Lacerta_multiple_lacvir1_roast.MAF >Lacerta_multiple_lacvir1_roast.SS; msa_view Lacerta_multiple_lacvir1_roast.SS.maf --in-format MAF --out-format SS > Lacerta_multiple_lacvir1_roast.SS; rm Lacerta_multiple_lacvir1_roast.SS.maf

 sed -e 2d Lacerta_multiple_lacvir1_roast.MAF >Lacerta_multiple_lacvir1_roast.MAF.edit

 phyloFit --msa-format MAF --tree "(anoCar2,(lacvir1,lacbil1)))" --out-root Lacerta_init Lacerta_multiple_lacvir1_roast.MAF.edit
# rm Lacerta_multiple_lacvir1_roast.MAF.edit 

printf "Created Lacerta_init model file" 

 msa_view --unordered-ss --aggregate lacvir1,lacbil1,anoCar2 --in-format SS --out-format SS new_codons_ss/*.codons.ss > Lacerta_all-cons.sites.ss; msa_view --unordered-ss --aggregate lacvir1,lacbil1,anoCar2,allMis1,galGal3 --in-format SS --out-format SS new_sites_ss/*.sites.ss > Lacerta_all-4d.sites.ss

printf "4dss files combined\n"


 phyloFit --tree "(anoCar2,(lacvir1,lacbil1))" --msa-format SS --out-root Lacerta_conserved-codon1 Lacerta_all-cons.sites.ss; phyloFit --tree "(anoCar2,(lacvir1,lacbil1))" --msa-format SS --out-root Lacerta_nonconserved-4d Lacerta_all-4d.sites.ss
 phyloBoot --read-mods Lacerta_conserved-codon1.0.mod,Lacerta_conserved-codon1.1.mod,Lacerta_conserved-codon1.2.mod,Lacerta_conserved-codon1.3.mod --output-average Lacerta_ave.cons.mod; phyloBoot --read-mods Lacerta_nonconserved-4d.0.mod,Lacerta_nonconserved-4d.1.mod,Lacerta_nonconserved-4d.2.mod,Lacerta_nonconserved-4d.3.mod --output-average Lacerta_ave.noncons.mod

printf "Cons and Noncons files created\n"


 phyloFit --tree "(anoCar2,(lacvir1,lacbil1))" --msa-format SS --out-root Lacerta_conserved-codon1 Lacerta_all-cons.sites.ss; phyloBoot --read-mods Lacerta_conserved-codon1.0.mod,Lacerta_conserved-codon1.1.mod,Lacerta_conserved-codon1.2.mod,Lacerta_conserved-codon1.3.mod --output-average Lacerta_ave.cons.mod; phyloBoot --read-mods Lacerta_nonconserved-4d.0.mod,Lacerta_nonconserved-4d.1.mod,Lacerta_nonconserved-4d.2.mod,Lacerta_nonconserved-4d.3.mod --output-average Lacerta_ave.noncons.mod

printf "Cons and Noncons files combined\n"


 for num in `seq 0.05 0.05 1`; do for num2 in `seq 5 3 50`; do printf NUM"\t"$num"\t"$num2"\n"; consEntropy $num $num2 Lacerta_ave.cons.mod Lacerta_ave.noncons.mod; done; done >consentropy_sq.txt
 for num2 in `seq 2 2 50`; do printf NUM"\t0.5\t"$num2"\n"; consEntropy 0.1 $num2 Lacerta_ave.cons.mod Lacerta_ave.noncons.mod; done >consentropy.txt_2sq
 printf "Check Entropy files\n\n"
 cat consentropy_sq.txt | grep -e NUM -e PIT | sed -e 's/.*=//' | tr "\n" " " | sed -e 's/ /\t/g' | sed -e 's/NUM/\nNUM/g' -e 's/\tbits/ bits/g' >consentropy.2.txt





ls new_mod_mafs/*.maf | parallel -j 30 sh run_new_phastcons_1.sh; ls new_Chunks/*.ss | parallel -j 30 sh run_new_phastcons_2.sh

ls new_Trees/*.cons.mod > Lacerta_cons.txt; ls new_Trees/*.noncons.mod > Lacerta_noncons.txt
for file in new_Trees/Lvir_*.cons.mod; do printf $file"\t"$(grep TREE $file | sed -e 's/.* (/(/' -e 's/,/ /g' -e 's/(+/ /g' -e 's/)+/ /g' | sed 's/ \s +/ /' | wc -w | awk '{print $0-1}')"\n"; done | sort -k2,2n >new_count_species

awk -F'\t' '{print >>"new_count_species."$2}' new_count_species; for file in new_count_species.[0-9]*; do awk -F'\t' '{print $1}' $file >new_phyloboot_$file; done
awk -F'\t' '{print >>"new_count_species."$2}' new_count_species; for file in new_count_species.[0-9]*; do awk -F'\t' '{print $1}' $file >new_phyloboot_$file; done
for file in new_phyloboot_*; do file2=$(echo $file | sed 's/.*\.//'); str=$(printf "'*"$file"'"); echo phyloBoot --read-mods $str --output-average Lacerta_ave.cons.mod_$file2 ; done; phyloBoot --read-mods '*new_phyloboot_new_count_species.1' --output-average Lacerta_ave.cons.mod_1; phyloBoot --read-mods '*new_phyloboot_new_count_species.2' --output-average Lacerta_ave.cons.mod_2 
for file in new_count_species.[0-9]*; do awk -F'\t' '{print $1}' $file | sed 's/\.cons\.mod/.noncons.mod/' >new_phyloboot_nc_$file; done; for file in new_phyloboot_nc_*.[1-2]; do file2=$(echo $file | sed 's/.*\.//'); str=$(printf "'*"$file"'"); cmd1=$(echo phyloBoot --read-mods $str --output-average Lacerta_ave.noncons.mod_$file2); eval $cmd1 ; done

ls new_Chunks/*.ss | parallel -j 30 sh run_new_phastcons_4.sh; cat new_Scores/Lvir_*.[1-2] >Lacerta_most-conserved.wig


cat new_Scores/Lvir_*.1 >Lacerta_most-conserved.1.wig; cat new_Scores/Lvir_*.2 >Lacerta_most-conserved.2.wig
for file in Lacerta_most-conserved.[1-2].wig; do /scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/wigToBigWig -clip $file vir.sizes $file.bigwig; done
/scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/bigWigMerge Lacerta_most-conserved.1.wig.bigwig Lacerta_most-conserved.2.wig.bigwig Lacerta_most-conserved.wig.bedgraph
/scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/bedGraphToBigWig Lacerta_most-conserved.wig.bedgraph vir.sizes Lacerta_most-conserved.bigwig
/scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/bigWigToWig Lacerta_most-conserved.bigwig Lacerta_most-conserved.wig

# /scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/bigWigToBedGraph Lacerta_most-conserved.4.wig.bigwig Lacerta_most-conserved.4.bedGraph


cat ./new_phylop_scores/Lvir_*.scores.wig >./Lacerta_phylop_all_scores.wig
/scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/wigToBigWig -clip Lacerta_phylop_all_scores.wig vir.sizes Lacerta_phylop_all_scores.bigwig
/scr/bloodymary/rohit/kentUtils-master/bin/linux.x86_64/bigWigMerge Lacerta_phylop_all_scores.bigwig Lacerta_phylop_all_scores.bedgraph


awk -F'\t' '$NF>0' Lacerta_phylop_base-scores.bed >Lacerta_phylop_base-scores.bed_positive
awk -F'\t' '$NF==0' Lacerta_phylop_base-scores.bed >Lacerta_phylop_base-scores.bed_neutral
awk -F'\t' '$NF<0' Lacerta_phylop_base-scores.bed >Lacerta_phylop_base-scores.bed_negative


awk -F'\t' '$NF<0' Lacerta_phylop_scores.txt  | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6,7 -o count,count_distinct,count_distinct,mean >./Lacerta_phylop_scores.txt_accl
awk -F'\t' '$NF==0' Lacerta_phylop_scores.txt  | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6,7 -o count,count_distinct,count_distinct,mean >./Lacerta_phylop_scores.txt_neut
awk -F'\t' '$NF>0' Lacerta_phylop_scores.txt  | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6,7 -o count,count_distinct,count_distinct,mean >./Lacerta_phylop_scores.txt_purf





