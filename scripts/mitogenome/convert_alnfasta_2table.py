#!/usr/bin/python3

### RUN as python convert_alnfasta_2table.py -i ALL.seb.aln -o ALL.seb.aln.list ###

from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--infile", "-i", type=str, required=True)
parser.add_argument("--outfile", "-o", type=str, required=True)
args = parser.parse_args()

def is_fasta(infile):
    with open(infile, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

is_fasta(args.infile)

outputfile=open(args.outfile, 'w')
print("sample","position","nucleotide", file=outputfile, sep='\t')

for seq_record in SeqIO.parse(args.infile, "fasta"):
    count=1
    for item in seq_record.seq:
        print(seq_record.id, count, item, file=outputfile, sep='\t')
        count+=1
