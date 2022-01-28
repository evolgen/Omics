#!/usr/bin/env python3

import sys
from Bio import AlignIO

fasta_seq_all = AlignIO.read(sys.argv[1], 'fasta')
fasta_seq_first = fasta_seq_all[0]

file_name1 = sys.argv[1]
orthogroup1 = file_name1[file_name1.find('OrthoGroup'):].split('.')[0]

if(len(list(fasta_seq_all)) != 1):
    print("Sequences = ", len(list(fasta_seq_all)), sep="")
    sys.exit (' Only one sequence input allowed (make it Human)')

align_length1 = int(len(fasta_seq_first))

print("Group", "Human_AA", "Aln_pos", "Human_pos", sep="\t")

actual_position = 0
for curr_base in range(0, align_length1):
    curr_base_exact = curr_base + 1
    list_curr_base_seq = []
    count_curr_base_seq = {}
    
    human_pos_exact = fasta_seq_first[curr_base]

    if human_pos_exact != "-":
        actual_position = actual_position + 1
    
    print(orthogroup1, human_pos_exact, curr_base_exact, actual_position, sep="\t")

