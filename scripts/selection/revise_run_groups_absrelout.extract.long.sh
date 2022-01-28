#!/usr/bin/bash

groupname=$1

all_branches_tree=$(fgrep -o -e 'Sebastes_' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ALN.cds.raw.tree.long | wc -l) ;

foreground_branches_tree=$(fgrep -o -e '{Foreground}' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ALN.cds.raw.tree.long | wc -l) ;

positive_branches_log=$(fgrep -e 'Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.log | sed -e 's/.* found \*\*//'  -e 's/\*\* branches under .*//') ;

foreground_branches_leg=$(fgrep -e 'Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.log | sed -e 's/.* among \*\*//'  -e 's/\*\* tested.*//') ;
#Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **1** branches under selection among **1** tested.

positive_species_names=$(fgrep -A 1000 -e 'Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.log | grep '^\* ' | sed -e '/^$/d' -e 's/^\* //' -e 's/, p-value =.*//' | tr "\n" "," | sed -e 's/,$/\n/') ;

positive_species_pvalues=$(fgrep -A 1000 -e 'Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.log | grep '^\* ' | sed -e '/^$/d' -e 's/.*, p-value =  //' | tr "\n" "," | sed -e 's/,$/\n/') ;

###printf "${groupname}\tAbsrel\tLong\t${all_branches_tree}\t${foreground_branches_tree}\t${positive_branches_log}\t${foreground_branches_leg}\t${positive_species_names}\t${positive_species_pvalues}\n" ;
#Group | Absrel | Long | All_species_tree | Foreground_species_tree | Positive_species_count | Foreground_species_count | Positive_species_names

echo $positive_species_names | sed -e 's/,/\n/g' >/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.speciesnames ;
echo $positive_species_pvalues | sed -e 's/,/\n/g' >/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.pvalues ;
cat /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.speciesnames | while read species; do fgrep -A 1000 '### Testing selected branches for selection' /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.log | fgrep -w "$species" | grep '^|'; done | sed -e 's/).*//' -e 's/.*|//' | tr -s " " | sed -e 's/^ //' -e 's/ /\t/' -e 's/%.*//' -e 's/( //' >/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.dnds

paste /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.dnds /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.pvalues /global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files/${groupname}/REVISE/ABSREL_raw/absrel_long.filt.speciesnames |
    awk -F'\t' -v groupname="$groupname" -v all_branches_tree="$all_branches_tree" -v foreground_branches_tree="$foreground_branches_tree" -v positive_branches_log="$positive_branches_log" -v foreground_branches_leg="$foreground_branches_leg" 'BEGIN{OFS="\t"} {print groupname,"Absrel","Long",all_branches_tree,foreground_branches_tree,positive_branches_log,foreground_branches_leg,$0}'

