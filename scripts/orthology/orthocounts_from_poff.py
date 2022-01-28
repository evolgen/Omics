#!/usr/bin/env python3

import sys
import csv
import argparse

#sebonly_tsv = "/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Geneloss/sebonly.lifted.poff.tsv"

if len(sys.argv) != 2:
    print("Please call one input file")
    sys.exit()

sebonly_tsv = sys.argv[1]

with open(sebonly_tsv) as tsv:
    speciesnames = tsv.readline().replace('.faa','').split("\t") # Read header, replace extension and split by tab
    speciesnames = [species1.replace('\n', '') for species1 in speciesnames]
    header_size = len(speciesnames)
    print ("GroupID", end='\t')
    count_header_elem1 = 0
    for curr_header1 in range(0, header_size):
        count_header_elem1 = count_header_elem1 + 1
        if (count_header_elem1 == header_size):
            print (speciesnames[curr_header1], end='')
        else:    
            print (speciesnames[curr_header1], end='\t')
    print()    
    count=0
    for ortholine1 in csv.reader(tsv, delimiter="\t"):
        count=count+1
        orthocount="Family{:06d}".format(count)
#            if(count > 5):
#                break
        number_of_species=len(ortholine1)
        print (orthocount, sep='\t', end='\t')
        curr_species1_elem = 0
        for curr_species1 in range(0, number_of_species):
            curr_species1_elem = curr_species1 + 1
            if (ortholine1[curr_species1] != "*" ):
                curr_species1_count = ortholine1[curr_species1].count(',')
                curr_species1_count = curr_species1_count + 1
            elif (ortholine1[curr_species1] == "*" ):
                curr_species1_count = 0
        
            if (curr_species1_elem != header_size):
                print (curr_species1_count, end='\t')
            else:
                print (curr_species1_count, end='\n')
            #print (speciesnames[curr_species1], curr_species1_count, sep='\t')    
    
print ("All counts done till ", orthocount, file=sys.stderr)

