#!/usr/bin/evn python3

import sys
import os
from collections import defaultdict

#====================================================================#

# read in the set of lines with dodgy genes on them in the embl file:

def read_lines_with_dodgy_genes(input_error_file):

    dodgy_gene_lines = set()
  
    fileObj = open(input_error_file, "r")
    for line in fileObj:
        line = line.rstrip()
        if line.startswith('ERROR:') and 'ERROR: Invalid AC *' not in line:
            # ERROR: Intron usually expected to be at least 10 nt long. Please check the accuracy. line: 699061-699062 of syphacia_muris.embl - syphacia_muris.embl
            # ERROR: Intron usually expected to be at least 10 nt long. Please check the accuracy. line: 699028 of syphacia_muris.embl - syphacia_muris.embl
            # ERROR: Abutting features cannot be adjacent between neighbouring exons. line: 1239561 of syphacia_muris.embl - syphacia_muris.embl
            temp = line.split("line: ")
            line = temp[1] # 699061-699062 of syphacia_muris.embl - syphacia_muris.embl 
            temp = line.split()
            line = temp[0] # 699061-699062
            temp = line.split('-')
            line = temp[0] # 699061
            dodgy_gene_lines.add(int(line))
    fileObj.close()

    return dodgy_gene_lines

#====================================================================#

# read in the set of lines with dodgy genes on them in the embl file:

def read_dodgy_genes(input_embl, dodgy_gene_lines):

    dodgy_genes = set()

    linecnt = 0  
    searching_for_gene_name = False
    gene_name = None
    fileObj = open(input_embl, "r")
    for line in fileObj:
        line = line.rstrip()
        linecnt += 1
        # FT   gene            complement(30891..35820)
        # FT                   /locus_tag="SMUV_LOCUS343"
        # FT                   /note="source:WormBase_imported"
        # FT                   /note="ID:gene:SMUV_0000034201"
        # FT                   /standard_name="SMUV_0000034201"
        if 'FT   gene' in line:
            searching_for_gene_name = True
            gene_name = None
        elif 'FT                   /standard_name="' in line:
            if searching_for_gene_name == True:
                assert(gene_name == None)
                temp = line.split('\"')
                gene_name = temp[1]
                searching_for_gene_name = False
        if linecnt in dodgy_gene_lines:
            # this is one of the lines that has a dodgy gene:
            assert(gene_name != None)
            dodgy_genes.add(gene_name)
    fileObj.close()

    return dodgy_genes

#====================================================================#

# write out the list of dodgy genes to the output file:

def write_output_file_of_dodgy_genes(dodgy_genes, output_file):

    outputfileObj = open(output_file, "w")
    for gene in dodgy_genes:
        output_line = "%s\n" % gene
        outputfileObj.write(output_line)
    outputfileObj.close()

    return

#====================================================================#

def main():
    
    # check the command-line arguments:
    if len(sys.argv) != 4 or os.path.exists(sys.argv[1]) == False or os.path.exists(sys.argv[2]) == False:
        print("Usage: %s input_embl input_error_file output_file" % sys.argv[0]) 
        sys.exit(1)
    input_embl = sys.argv[1] # e.g. syphacia_muris.embl
    input_error_file = sys.argv[2] # e.g. syphacia_muris_VAL_ERROR.txt  
    output_file = sys.argv[3] # output file for putting a list of dodgy genes

    # read in the set of lines with dodgy genes on them in the embl file:
    dodgy_gene_lines = read_lines_with_dodgy_genes(input_error_file)

    # read in the embl file to find the genes highlighted as dodgy by the ENA's EMBL file validator:
    dodgy_genes = read_dodgy_genes(input_embl, dodgy_gene_lines)

    # write out the list of dodgy genes to the output file:
    write_output_file_of_dodgy_genes(dodgy_genes, output_file)

    print("FINISHED\n")

#====================================================================#

if __name__=="__main__":
    main()

#====================================================================#
