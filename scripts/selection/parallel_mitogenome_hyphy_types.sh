#!/usr/bin/bash

set -e -o pipefail

module load trimal mafft raxml gcc muscle java macse;
#conda activate hyphy;

#for type in All Internal Leaves; do 

#namen="aleutianus";
type="All";

branch="${type}";
out_dir="$PWD/hyphy_absrel_${type}";
tree_file="/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.newick";
#"/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/species_tree.aleutianus.newick";
sequence_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/sequence_files";

mkdir -p ${out_dir};

genefile=$@

        printf "\t${genefile}\n";
        gene_name=$(echo $genefile | sed  -e 's/.FASTA$//' -e 's/.FAS$//' -e 's/.FA$//' -e 's/.fasta$//' -e 's/.fas$//' -e 's/.fa$//');

        cat /global/scratch2/rohitkolora/Rockfish/Genomes/Selection/mitogenome/list_genomes_seb.txt | 
            sed -e 's/$/\/mitogenome\/MITOS2\/CDS\//' |
            while read genomedir; do
                speciesname=$(echo $genomedir | sed -e 's/.*\/output\///' -e 's/\/.*//');
            done

        mkdir -p ${out_dir}/${gene_name}_absrel ${out_dir}/${gene_name}_busted ; 
        mkdir -p ${out_dir}/${gene_name}_fel ${out_dir}/${gene_name}_relax ${out_dir}/${gene_name}_meme ;
##        cd ${out_dir}/${gene_name}_absrel;
##        hyphy aBSREL --code Vertebrate-mtDNA --alignment ${sequence_dir}/${genefile}.aln.nex --tree ${tree_file} --branches ${branch} --output ${genefile}.json --save-fit ${genefile}.fit 1>${genefile}.log 2>&1
##        cd ${out_dir}/${gene_name}_busted;
##        hyphy busted --code Vertebrate-mtDNA --alignment ${sequence_dir}/${genefile}.aln.nex --tree ${tree_file} --branches ${branch} --output ${genefile}.json --save-fit ${genefile}.fit 1>${genefile}.log 2>&1
##        cd ${out_dir}/${gene_name}_fel;
##        hyphy fel --code Vertebrate-mtDNA --alignment ${sequence_dir}/${genefile}.aln.nex --tree ${tree_file} --branches ${branch} --output ${genefile}.json --save-fit ${genefile}.fit 1>${genefile}.log 2>&1
##        #cd ${out_dir}/${gene_name}_relax;
##        #hyphy busted --code Vertebrate-mtDNA --alignment ${sequence_dir}/${genefile}.aln.nex --tree ${tree_file} --branches ${branch} --output ${genefile}.json --save-fit ${genefile}.fit 1>${genefile}.log 2>&1
        cd ${out_dir}/${gene_name}_meme;
        hyphy meme --code Vertebrate-mtDNA --alignment ${sequence_dir}/${genefile}.aln.nex --tree ${tree_file} --branches ${branch} --output ${genefile}.json --save-fit ${genefile}.fit 1>${genefile}.log 2>&1


