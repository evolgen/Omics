#!/usr/bin/bash

set -e -o pipefail

query="/global/scratch2/rohitkolora/Other/Hummingbird/broadtail.fasta"
reference="/global/scratch2/rohitkolora/Other/Hummingbird/Anna/GCF_003957555.fasta"

module load minimap2 samtools bedtools ;
printf "\tRunning minimap - First\n";

minimap2 -t $(eval nproc) -x asm10 --score-N 2 -g 1000 ${reference} ${query} | cut -f1-12 >paf1 ; 
bash ~/RGP/scripts/rearrangements/pairsums_paf2.sh paf1 paf1.sum paf1.top ; 

printf "\tExtracted top sequence pairs\n" ;
grep '^>' ${reference} | sed -e 's/^>//' -e 's/ .*//' |
  while read name; do
  grep -w -e "$name" paf1.top | cut -f1 ;
      done | awk '!x[$1]++' |
      xargs samtools faidx ${query} > querygenlisted1.fasta ;                            # GET the query genome - first run
  printf "\tEditing fasta headers\n" ;
  bash ~/RGP/scripts/edit_chrheader.sh ${reference} refgenome ;
  bash ~/RGP/scripts/edit_chrheader.sh querygenlisted1.fasta qrygenome;
  echo "###REF-names###" >syri1.sp1.sp2.chrnames;
  paste <(grep '^>' ${reference}) <(grep '^>' refgenome) >>syri1.sp1.sp2.chrnames;
  echo "###QRY-names###" >>syri1.sp1.sp2.chrnames;
  paste <(grep '^>' querygenlisted1.fasta) <(grep '^>' qrygenome) >>syri1.sp1.sp2.chrnames;       # CREATE pair aliases
  printf "\tExecuting minimap for 1st syri iteration - Second\n";
  minimap2 -t $(eval nproc) -ax asm10 --score-N 2 -g 1000 --cs=long -Lc refgenome qrygenome >align.sp1.sp2.sam; 
  printf "\tConverting it to bam\n";        
  bash ~/RGP/scripts/sam2bam_withmatch.sh align.sp1.sp2.sam align.sp1.sp2.bam ; 
  rm -f align.sp1.sp2.sam ;     
  printf "\tExecuting SYRI - Initial\n";
  printf "#!/usr/bin/bash\n\n" >iter1.syri.run.bash ;
  echo "python3 $PATH_TO_SYRI -c align.sp1.sp2.bam -r refgenome -q qrygenome -k -F B --no-chrmatch --log WARN --novcf --nosr --nosv --nosnp --nc $(eval nproc) --prefix iter1.syri.sp1.sp2. 1>log_syri1 2>&1" >>iter1.syri.run.bash ;
  bash iter1.syri.run.bash;    
  samtools faidx qrygenome ;    
  grep -e 'has high fraction of inverted alignments' log_syri1 |
      sed -e 's/.* query genome (//' -e 's/)\. Ensure that.*//' >log_syri1_rcnames;    # GET probable Rev-Compl
  printf "" >qrygenome.listedby.sp1 ;
  printf "\tReordering the query genome\n";
  bash ~/RGP/scripts/rearrangements/convert_rc_4mlog-qrygen.sh log_syri1_rcnames qrygenome qrygenome.listedby.sp1 ; 
  samtools faidx qrygenome.listedby.sp1 ;
  cut -f1 qrygenome.listedby.sp1.fai >qrygenome.listedby.sp1.names ;               # GET names only for RC
  cut -f1 qrygenome.fai |
      grep -w -v -f qrygenome.listedby.sp1.names |
      xargs samtools faidx qrygenome >>qrygenome.listedby.sp1 ;                      # ADD rest of the sequences
  samtools faidx qrygenome.listedby.sp1 ;
  cut -f1 qrygenome.fai |
      xargs samtools faidx qrygenome.listedby.sp1 >qrygenome.listedby.sp1.edit ; 
  samtools faidx qrygenome.listedby.sp1.edit ;
  printf "\tCreating intermediate files before SYRI\n";
  bash ~/RGP/scripts/edit_chrheader.sh qrygenome.listedby.sp1.edit qrygenome_fin ;
  echo "###REF-names###" >syri2.sp1.sp2.chrnames;
  paste <(grep '^>' ${reference}) <(grep '^>' refgenome) >>syri2.sp1.sp2.chrnames;
  echo "###QRY-names###" >>syri2.sp1.sp2.chrnames;
  paste <(grep '^>' ${query}) <(grep '^>' qrygenome_fin) >>syri2.sp1.sp2.chrnames;  # CREATE PAIR aliases
  printf "\tExecuting minimap for SYRI - THIRD \n";
  minimap2 -t $(eval nproc) -ax asm10 --score-N 2 -g 1000 --cs=long -Lc refgenome qrygenome_fin >align.sp1.sp2.sam;
  printf "\tBAMing for SYRI - 2nd iteration\n";        
  bash ~/RGP/scripts/sam2bam_withmatch.sh align.sp1.sp2.sam bam ; 
  rm -f align.sp1.sp2.sam ;
  printf "\tExecute SYRI - Part II \n" ;

  module load minimap2 samtools bedtools ;
  printf "\tExecuting SYRI - FINALLY\n" ;
  python3 $PATH_TO_SYRI -c bam -r refgenome -q qrygenome_fin -b 5 -k -F B --no-chrmatch --nosnp --lf Final.syri.LOG --nc $(eval nproc) --cigar --prefix Final.sp1.sp2. 1>LOG_syri2 2>&1;
  bash ~/RGP/scripts/rearrangements/execute_synt-replacer.sh Final.sp1.sp2.syri.out syri1.sp1.sp2.chrnames Final.sp1.sp2.syri.synt.list ;
  bash ~/RGP/scripts/rearrangements/execute_syriout-namereplacer.sh Final.sp1.sp2.syri.out syri1.sp1.sp2.chrnames Final.sp1.sp2.syri.out.TSV ;

