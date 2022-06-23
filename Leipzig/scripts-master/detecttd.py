#! /usr/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
# Name: detecttd.py
# Purpose:Detecting tandem duplications in NGS short reads
#
# Author: Fabian Grandke
#
# Created: 2012/10/25
# Copyright: (c) 2012
# Licence: <your licence>
# Homepage: http://sourceforge.net/projects/detecttd/
#----------------------------------------------------------------------------- 


import argparse
import os
import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import operator

### Try to import time module (to create temporary file names)
try:
  from time import time
  time_imp = True
except ImportError:
  time_imp = False 
  
def main():
  #Command Line Parser:
  parser = argparse.ArgumentParser(
	    description='Program to find Tandem Duplications in NGS reads' )
  parser.add_argument('-i','--input', metavar = 'FILE', 
	    help = 'Inputfile FASTA', required=True)
  parser.add_argument('-o', '--output', metavar = 'FILE', 
	    help = 'Prefix for outputfiles', required=True)
  #parser.add_argument('-f', '--format', metavar = 'STRING',
	    #choices = ['sanger','solexa','illumina'], default = 'sanger',
			#help = 'Type of FASTQ format ("sanger", "solexa" or "illumina")')
  parser.add_argument('-m', '--multiplex', default = 0, metavar = 'INT', 
	    help = 'Number of nucleotides to cut off, in case of multiplexing')
  parser.add_argument('-t', '--temp', default = "/tmp", 
	    help = 'Location to store temporary files', metavar = 'DIR')
  parser.add_argument('-s', '--seed', help = 'Length of seed in BLAST',
	    type = int, metavar='INT', default = 15)
  parser.add_argument('-d', '--identity', help = 'Minimal identity of blast hits',
	    type = int, metavar='INT', default = 95)
  parser.add_argument('-r', '--reference', help = 'If given, the modified reads will be blasted against this reference.',
	    metavar='FILE', default = '')	    
  #parser.add_argument('-l', '--lonely', help = 'Map duplicated region without surroundings.',
	    #metavar='INT', default = 0)	    
  parser.add_argument('-p', '--penalty', help = 'Penalty for blast command.',
	    metavar='INT', default = -4)		    
  args = parser.parse_args()
  
  if(int(args.penalty) > 0):
    sys.exit('penalty must be <= 0')
  ###Check Inputfile
  fileName, fileExtension = os.path.splitext(args.input)
  if fileExtension == '.fasta' or fileExtension == '.fa':
    format = 'fasta'
  elif fileExtension == '.fastq' or fileExtension == '.fq':
    format = 'fastq' + "-" + args.format
  else:
    sys.exit('Inputfile is not FASTA or FASTQ')
    
  ### Temporary files
  if time_imp:
    query = args.temp + "/" + os.path.splitext(os.path.basename(args.input))[0] + str(time()) + ".seq"
    blast = args.temp + "/" + os.path.splitext(os.path.basename(args.input))[0] + str(time()) + ".blast"
  else:
    query = args.temp + "/" + os.path.splitext(os.path.basename(args.input))[0] +  ".seq"
    blast = args.temp + "/" + os.path.splitext(os.path.basename(args.input))[0] +  ".blast"    
  
  ### Output files
  orig = args.output + "_orig." + format.split('-')[0]
  smap = args.output + "_smap." + format.split('-')[0]
  out = args.output + ".tdf"
  result = args.output + ".out"
  origh = open(orig, "w")
  smaph = open(smap, "w")
  outh = open(out, "w")
  ###Print header to tdf
  outh.write("###### Output of Tandem Duplication Finder ###\n" +
	      "### Read ID: ID of the TD-containing read\n" +
	      "### Length: Number of duplicated nucleotides\n" +
	      "### Gapsize: Number of nucleotides between duplications\n" +
	      "### Identity: Percentage of identity between duplications\n" +
	      "### Frameshift: Shift of reading frame produced by duplication and gap (0, 1 or 2)\n" +
	      "### Read ID\tLength\tGapsize\tEdge\tIdentity\tFrameshift\n" )
  ###Open input
  input = open(args.input, "rU")
  rec_counter = 0
  td_counter = 0
  td_length = []
  success = 0
  edges = [0,0,0,0]
  reads_in_orf = 0
  mapcount = 0
  pos_list = []
  for record in SeqIO.parse(input, format) :
    rec_counter = rec_counter + 1
    queryh = open(query, "w")
    SeqIO.write(record, queryh, "fasta")
    queryh.close()
    ### Run BLAST
    blastn_command = NcbiblastnCommandline(query = query ,
    subject = query , outfmt = 5, word_size = args.seed, 
		perc_identity = args.identity ,out = blast,
		penalty = args.penalty)
    stdout, stderr = blastn_command()
    ### Parse blast output
    blasth = open(blast)
    blast_records = NCBIXML.parse(blasth)
    blast_record = blast_records.next()
    ### Try and except out of IndexError
    try:
      num_ali = blast_record.descriptions[0].num_alignments
    except IndexError:
      continue
    ### Selfhit + tandem duplication found
    if num_ali > 2:
      for alignment in blast_record.alignments:
	### new_align = True
	found = False
	for hidx, hsp in enumerate(alignment.hsps):
	  for hidx2, hsp2 in enumerate(alignment.hsps):
	    if hidx2 <= hidx:
	      continue
	    if (hsp.query_start == hsp2.sbjct_start and 
		    hsp2.query_start == hsp.sbjct_start and 
		    hsp.strand[0] == hsp.strand[1]):			#both hsps on the same strand
	      found = True
	      ###Overlap detection, start and end of td are equal and end will be shifted
	      if(hsp.sbjct_end >= hsp2.sbjct_start):
		diff = hsp.sbjct_end - hsp2.sbjct_start + 1
		hsp.query_end -= diff
		hsp2.query_end -= diff
		hsp.query = hsp.query[:( -diff)]
		hsp2.query = hsp2.query[:( -diff)]
		hsp.sbjct_end -= diff
		hsp2.sbjct_end -= diff
		hsp.sbjct = hsp.sbjct[:( -diff)]		
		hsp2.sbjct = hsp2.sbjct[:( -diff)]
		hsp.match = hsp.match[:( -diff)]
		hsp2.match = hsp2.match[:( -diff)]
		
	      pos_list.append(sorted([hsp.query_start, hsp.query_end, hsp2.query_start, hsp2.query_end]))
	      td_counter = td_counter + 1	      
	      ### indicate tandem duplications
	      td_marker = ''.join(' ' for i in range(1,hsp.sbjct_start)) + '>'
	      td_marker += ''.join('-' for i in range(1,len(hsp.sbjct)-1)) + '<'
	      td_marker += ''.join(' ' for i in range(hsp.sbjct_end,hsp2.sbjct_start-1)) + '>'
	      td_marker += ''.join('-' for i in range(1,len(hsp2.sbjct)-1)) + '<'
	      ### print sequence with tandem duplications in lower case
	      start = str(record.seq[0:hsp.sbjct_start-1])				#start
	      f_td_orig = str(record.seq[hsp.sbjct_start-1:hsp.sbjct_end])
	      #print f_td_orig
	      f_td_mod = hsp.sbjct
	      #print f_td_mod
	      gap = str(record.seq[hsp.sbjct_end : hsp2.sbjct_start-1])		#gap
	      s_td_orig = str(record.seq[hsp2.sbjct_start-1:hsp2.sbjct_end])
	      s_td_mod = hsp2.sbjct
	      end = str(record.seq[hsp2.sbjct_end:])				#end
	      case_sen = start.upper() + f_td_mod.lower() + gap.upper() + s_td_mod.lower() + end.upper()
	      ### build up string indicating (mis-)matches and gap
	      matcher = ''.join(' ' for i in range(1,hsp.sbjct_start)) + hsp.match
	      matcher += ''.join(' ' for i in range(hsp.sbjct_end,hsp2.sbjct_start-1)) + hsp2.match
	      ### build up string showing alignment
	      td_seqs = ''.join(' ' for i in range(1,hsp.sbjct_start)) + hsp.query.lower()
	      td_seqs += ''.join(' ' for i in range(hsp.sbjct_end,hsp2.sbjct_start-1)) + hsp2.query.lower()

	      ### Output original sequence in original format
	      SeqIO.write(record, origh, format)
	      ### Print output of to output file
	      identity_val = float(int((float(min(hsp.identities, hsp2.identities)) / float(hsp.align_length))*10000))/100.0	      
	      outh.write(record.id + "\t" + str(len(f_td_mod)) + "\t" + str(len(gap)) + "\t" + 
		    str(identity_val) + "\t" + str((len(s_td_mod) + len(gap)) % 3) + "\n+ " +
		    td_marker + "\n+ " + case_sen + "\n+ " + matcher+ "\n+ " + td_seqs + "\n")		    
	      break
	  if found: break

    #########################################
  input.close()
  origh.close()

  if(args.reference != ''):
    ##############BLAST ORIGINAL SEQUENCES and expect to sbjct-overlapping, query-disjunct hits ->megablast
    blastn_command = NcbiblastnCommandline(query = orig,
      subject = args.reference , outfmt = 5, word_size = 11, 
	      out = blast, perc_identity = 0.95, )
    stdout, stderr = blastn_command()
    blasth = open(blast)
    blast_records = NCBIXML.parse(blasth)
    td_list = [] #seq, start, end, counter
    for index1, blast_record in enumerate(blast_records):  
      try:
	num_ali = blast_record.descriptions[0].num_alignments
      except IndexError:
	continue
      if(num_ali >= 1): 	
        mapcount += 1
        #TODO
        smaph.write(blast_record.query + "\n")
	alignment = blast_record.alignments[0]
	hsp1 = alignment.hsps[0]	
	sbjct_pos = [hsp1.sbjct_start, hsp1.sbjct_end]
	  	  
	hsp1_sbjct_strand = 1
	if(hsp1.sbjct_start > hsp1.sbjct_end):
	  sbjct_pos[0] = hsp1.sbjct_start * (-1)
	  sbjct_pos[1] = hsp1.sbjct_end * (-1)
	if(hsp1.query_start <= pos_list[index1][0] and hsp1.query_end >= pos_list[index1][1]): #first td is located within hsp1 and thus, could be mapped
	  abs_start = pos_list[index1][0]
	  start_pos = pos_list[index1][0] - hsp1.query_start					#start position in hsp1.query = first letter of first td
	  abs_end = min(hsp1.query_end,pos_list[index1][2]-1)
	  end_pos = abs_end - hsp1.query_start							#end position in hsp1.query
	  gap = pos_list[index1][2] - abs_end -1						#gap between end of first copy and second copy start
	elif(hsp1.query_start <= pos_list[index1][2] and hsp1.query_end >=pos_list[index1][3]):  #second td is located within hsp1 and thus, could be mapped
	  abs_start = max(hsp1.query_start, pos_list[index1][1]+1)
	  start_pos = abs_start - hsp1.query_start						#start position in hsp1.query	  
	  abs_end = pos_list[index1][3]
	  end_pos = pos_list[index1][3] - abs_start + start_pos				#end position in hsp1.query
	  gap = abs_start - pos_list[index1][1] -1						#gap between end of first copy and second copy start[2]) + "\t" + str(pos_list[index1][3])
	else:
	  mapcount -= 1
	  continue
	dupl_seq = str(hsp1.query)[start_pos:end_pos + 1]
	sbjct_pos[0] += start_pos
	sbjct_pos[1] -= end_pos
	td_length.append(len(dupl_seq))
	
	frame_shift = (len(str(dupl_seq)) + gap) % 3
	if frame_shift == 0:
	  reads_in_orf += 1
	td_list.append([1, sbjct_pos[0], sbjct_pos[1], gap, dupl_seq, frame_shift])
    ###Unique list and count duplications
    indices = range(len(td_list))
    for index, td in enumerate(td_list):
      uni_list = td_list
      for index2, td2 in enumerate(td_list):
	if(index2 <= index):
	  continue
	if(td_list[index][4] == td2[4]):
	  indices[index2] = -1
	  td_list[index][0] += 1 
    uni_list = []
    for indi in indices:
      if indi >=0 :
	uni_list.append(td_list[indi])
    uni_list.sort(key=operator.itemgetter(0,5))
    ###Print data to outputfile
    resulth = open(result, "w")
    resulth.write("#Reads\tLength\tStart\tGapsize\tFrame\tSequence\n")
    for index3, td3 in enumerate(uni_list):
      resulth.write(str(td3[0]) + "\t" + str(len(str(td3[4]))) + "\t" + str(td3[1]) + "\t" +
       str(td3[3]) + "\t" + str(td3[5]) + "\t" + str(td3[4]) + "\n")

    resulth.close()
  ### Write summary to output file
  try:			#if no td's have been found td_length is empty and throws errors, when accessed by max()
    max(td_length)
  except ValueError:
    td_length.append(0)
  
  edge_string = str(edges[0]) + " / " + str(edges[1]) + " / " + str(edges[2]) + " / " + str(edges[3]) 
  avg_length = str(float(sum(td_length)) / len(td_length))
  outh.write("# Sequences inspected: " + str(rec_counter) + "\n# Duplication candidates detected: " + str(td_counter) + 
  "\n# Max TD Length: " + str(max(td_length)) + "\n# Min TD Length: " + str(min(td_length)) + "\n# Average TD Length: " + 
  avg_length + "\n# Mapped reads: " + str(mapcount)+ "\n# TD's in readingframe: " + str(reads_in_orf) + "\n")
  outh.close()  
  os.system("rm -f " + query + " " + blast)  

  ##########################################

  
if __name__ == "__main__":
    main()