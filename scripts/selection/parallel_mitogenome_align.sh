#!/usr/bin/bash

set -e -o pipefail

module load trimal mafft raxml gcc muscle java macse;
#conda activate hyphy;

mkdir -p $PWD/sequence_files $PWD/hyphy_busted;
workingdir="/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome"

genefile=$@

#cat /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/list_mitogenes.txt | while read genefile; do
        printf "" >$PWD/sequence_files/${genefile};
        printf "${genefile}\t:\t";

        gene_name=$(echo $genefile | sed  -e 's/.FASTA$//' -e 's/.FAS$//' -e 's/.FA$//' -e 's/.fasta$//' -e 's/.fas$//' -e 's/.fa$//');
        cat /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/list_genomes_sebtree.txt | 
            sed -e 's/$/\/mitogenome\/MITOS2\/CDS\//' |
            while read genomedir; do
                speciesname=$(echo $genomedir | sed -e 's/.*\/output\///' -e 's/\/.*//');
                printf "${speciesname}\t";
                printf ">${speciesname}\n" >>$PWD/sequence_files/${genefile};
                echo ${genomedir}/${genefile}
                sequence_string=$(cat ${genomedir}/${genefile} | grep -v '^>' | fold -w3 | sed -e '/^$/d' -e '$s/TAA$/NNN/' -e '$s/TAG$/NNN/' -e '$s/TGA$/NNN/' -e '$s/taa$/NNN/' -e '$s/tag$/NNN/' -e '$s/tga$/NNN/' | awk 'length($0)==3' | tr "\n" " " | sed -e 's/ //g' -e 's/$/\n/');
                sequence_length=$(printf ${sequence_string} | awk '{ print length($0)%3 }');

                if [ "${sequence_length}" -eq "0" ]; then
                    echo ${sequence_string} >>$PWD/sequence_files/${genefile};
                fi
                if [ "${sequence_length}" -eq "1" ]; then
                    echo ${sequence_string} | sed -e 's/[a-zA-Z]$//' >>$PWD/sequence_files/${genefile};
                fi
                if [ "${sequence_length}" -eq "2" ]; then
                    echo ${sequence_string} | sed -e 's/[a-zA-Z][a-zA-Z]$//' >>$PWD/sequence_files/${genefile};
                fi    
            done
            echo;
        muscle -diags -in ${workingdir}/sequence_files/${genefile} -out ${workingdir}/sequence_files/${genefile}.aln;
        trimal -gt 1 -in ${workingdir}/sequence_files/${genefile}.aln -nexus -out ${workingdir}/sequence_files/${genefile}.nex; 

#        printf '' >$PWD/hyphy_busted/${genefile}.bf;
#        printf 'inputRedirect = {};\n' >>$PWD/hyphy_busted/${genefile}.bf;
#        printf 'inputRedirect["01"] = "Vertebrate mtDNA";\n' >>$PWD/hyphy_busted/${genefile}.bf;
#        printf 'inputRedirect["02"] = "sequence_files/${genefile}.nex";\n' >>$PWD/hyphy_busted/${genefile}.bf;
#        printf 'inputRedirect["03"] = "/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick";\n' >>$PWD/hyphy_busted/${genefile}.bf;
#        printf 'inputRedirect["04"] = "All";\n' >>$PWD/hyphy_busted/${genefile}.bf;
#        printf 'inputRedirect["05"] = "$PWD/hyphy_busted/${genefile}.result";\n' >>$PWD/hyphy_busted/${genefile}.bf;
#        hyphy ${genefile}.bf 1>${genefile}.log 2>&1;

        mkdir -p ${workingdir}/hyphy_busted/${gene_name};
        cd ${workingdir}/hyphy_busted/${gene_name};
        hyphy busted --code Vertebrate-mtDNA --alignment ${workingdir}/sequence_files/${genefile}.nex --tree /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick --branches All --output ${genefile}.json --save-fit ${genefile}.fit 1>${genefile}.log 2>&1


#done

## 10 4 1 file Y 1 d ## Hyphy options
#inputRedirect = {};
#inputRedirect["01"] = "Universal";
#inputRedirect["02"] = "/path/to/sequence/file";
#inputRedirect["03"] = "/path/to/tree/file";
#inputRedirect["04"] = "All";
#inputRedirect["05"] = "";

