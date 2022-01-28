#!/usr/bin/env python3

import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

alignment = AlignIO.read(sys.argv[1], 'fasta')
file_name1 = sys.argv[1]
orthogroup1 = file_name1[file_name1.find('OrthoGroup'):].split('.')[0]

summary_align = AlignInfo.SummaryInfo(alignment)

#consensus_seq = summary_align.dumb_consensus(threshold=float(sys.argv[2]), ambiguous='-')
consensus_seq = summary_align.gap_consensus(threshold=float(sys.argv[2]), ambiguous='-')

consensus_seq = list(consensus_seq.strip())

align_length1 = int(len(consensus_seq))

print("Group", "Aln_pos", "backspecies_consensus_aa",sep="\t")
for curr_base in range(0, align_length1):
    curr_base = int(curr_base)
    curr_base_exact = curr_base + 1
    print(orthogroup1, curr_base_exact, consensus_seq[curr_base], sep="\t")	

