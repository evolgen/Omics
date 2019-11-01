#!/usr/bin/env python3
import regex as re
import pdb
import numpy as np
import sys
import pandas as pd

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

###l = open("med_chimp2human.paf").readline()
for l in open(sys.argv[1],'r').readlines():

    sl = l.split("\t")
    ##sld = pd.DataFrame(sl) 
    ##slp = (sld[X11] == "60")

    if sl[11] != '60':    # Mapping quality 60 filter
        continue

    c = sl[-2]
    c = c.split(":")[-1]

    #cs = [int(i) for i in re.compile("M|D|I").split(c)]
    cs = re.compile("M|D|I").split(c)
    cs = np.array([int(i) for i in cs if i!=""])
    strs = np.array([s for s in re.compile("[0-9]*").split(c) if s != ""])

    #print("paf col 11 = %s"%sl[10])
    #print("sum of bases= %s"%np.sum(cs))
    
    aligned = np.sum(cs[strs=="M"]) #matches and mismatches
    matches = int(sl[9])    #total match bp
    block_len = int(sl[10]) #length of alignment block
    mismatches= aligned - matches   #total mismatch bp
    n_gap = np.sum(strs=="I")+np.sum(strs=="D") #Number gaps due to Insertions and deletions
    gap_len = np.sum(cs[strs=="I"])+np.sum(cs[strs=="D"])   #Total bases in Insertions and deletions

    #print(mismatches / (block_len - gap_len))
    #print((mismatches+n_gap) / (block_len - gap_len))
    heng_ID = (matches / (block_len - gap_len+n_gap))  #Heng Li's gap compressed divergence
    heng_Div = 1-heng_ID

    print(str(sl[0]) + "\t" + str(sl[2]) + "\t" + str(sl[3]) + "\t" + str(sl[5]) + "\t" + str(sl[7]) + "\t" + str(sl[8]) + "\t" + "\t" + str(aligned) + "\t" + str(matches) + "\t" + str(mismatches) + "\t" + str(heng_Div) + "\t" + str(mismatches / aligned))

#pdb.set_trace()
    
