#!/usr/bin/bash

#set -e

# conda activate py2.7
for file1 in SraExperimentPackage\ \(*.xml; do 
  
  ID=$(fgrep -m1 accession "$file1" | sed -e 's/.*accession="SRX/SRX/' -e 's/".*//' -e 's/^<EXPERIMENT_PACKAGE>.*"//') ; 
  printf "Processing started for $file1\t$ID\n";
  mkdir -p $ID; cp "$file1" ./${ID}/SraExperimentPackage.xml; 
  cd ${ID}; 
  python2 ~/code/seqlib/seqlib/sra/generate_download_config.py --fn_xml=SraExperimentPackage.xml
  cd ..;
  printf "$file1:\tDONE\n"

done


