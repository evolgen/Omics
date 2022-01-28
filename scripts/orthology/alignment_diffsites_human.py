#!/usr/bin/env python3

import sys
from Bio import SeqIO
import decimal
import collections

if len(sys.argv) != 2:
    print("Please call one input file")
    sys.exit()

record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
file_name1 = sys.argv[1]
orthogroup1 = file_name1[file_name1.find('OrthoGroup'):].split('.')[0]

#record_dict = SeqIO.to_dict(SeqIO.parse('/global/scratch2/rohitkolora/Rockfish/Genomes/Provean/long_PSG_dnarepairGO/sorted.OrthoGroup1116.FAA', 'fasta'))
#file_name1 = "/global/scratch2/rohitkolora/Rockfish/Genomes/Provean/long_PSG_dnarepairGO/sorted.OrthoGroup1116.FAA"
#orthogroup1 = file_name1[file_name1.find('OrthoGroup'):].split('.')[0]

all_species =  list([*record_dict])

seq_lengths1 = [len(val1) for val1 in record_dict.values()]
uniq_seq_lengths = sorted(set(seq_lengths1))
size_uniq_seq_lengths = len(uniq_seq_lengths)
if(size_uniq_seq_lengths != 1):
    sys.exit (' The sequences must be an MSA of equal length')

if(len(all_species) != 2):
    sys.exit (' We need atleast two-sequences MSA first')
    
first_species = list(record_dict.keys())[0]
align_length1 = len(record_dict[first_species])

#print(first_species, align_length1)

curr_seq_dict = dict([(x1,record_dict[x1]) for x1 in list(record_dict.keys())])

seqnames1 = "\t".join(list(curr_seq_dict.keys()))

human_description = list(curr_seq_dict.items())[0][-1].description
human_description = human_description[human_description.find(' GN='):]
gene_name = human_description.split(" ")[1].replace("GN=","")

print("#Gene", "OrthoGroup", "Human_seq", "Rockfish_seq", 
      "Human_AA", "Rockfish_AA", "Aln_pos", "Human_pos", 
      "Provean_input", sep="\t")

actual_position = 0
for curr_base in range(0, align_length1):
    curr_base_exact = curr_base + 1
    list_curr_base_seq = []
    count_curr_base_seq = {}
    for key1, value1 in curr_seq_dict.items():
        list_curr_base_seq.append(curr_seq_dict[key1][curr_base])
    human_pos_exact = list(list_curr_base_seq[0])    
    if "-" not in human_pos_exact:
        actual_position = actual_position + 1     
    if "-" in list(list_curr_base_seq):
        continue
    if len(sorted(set(list_curr_base_seq))) == 1:
        continue
    provean_currbase = str(list_curr_base_seq[0]) + str(actual_position) + str(list_curr_base_seq[1])    
    print(gene_name, orthogroup1, seqnames1, list_curr_base_seq[0], list_curr_base_seq[1], 
          curr_base_exact, actual_position, 
          provean_currbase, sep="\t")    

