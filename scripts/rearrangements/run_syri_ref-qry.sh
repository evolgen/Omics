#!/usr/bin/bash

set -e -o pipefail

# bash ~/RGP/scripts/rearrangements/run_syri_ref-qry.sh ref.fasta query.fasta outdir

reference=$1
query=$2
fulldir=$3

mkdir -p $fulldir && cd $fulldir ;

module load minimap2 samtools bedtools ;
printf "\tRunning minimap - First\n";

if [[ ! -f "paf1.top" ]]; then
    minimap2 -t $(eval nproc) -x asm10 --score-N 2 -g 1000 ${reference} ${query} | cut -f1-12 >paf1 ; 
    bash ~/RGP/scripts/rearrangements/pairsums_paf2.sh paf1 paf1.sum paf1.top ; 
fi

if [[ ! -e "${reference}.fai" ]]; then
    samtools faidx ${reference} ;
fi
if [[ ! -e "${query}.fai" ]]; then
    samtools faidx ${query} ;
fi    

printf "\tExtracted top sequence pairs\n" ;

printf "\tCreating first bam\n" ;
if [[ ! -f "align.sp1.sp2.bam" ]]; then
  cut -f1 ${reference}.fai | 
    fgrep -w -f - paf1.top |
    awk -F'\t' '!x[$1]++ {print $1}' >qry.list.names ;
      
  cat qry.list.names | xargs samtools faidx ${query} > querygenlisted1.fasta ;                            # GET the query genome - first run
  
  printf "\tEditing fasta headers\n" ;
  bash ~/RGP/scripts/edit_chrheader.sh ${reference} refgenome ;
  bash ~/RGP/scripts/edit_chrheader.sh querygenlisted1.fasta qrygenome;
  echo "###REF-names###" >syri1.sp1.sp2.chrnames;
  paste <(grep '^>' ${reference}) <(grep '^>' refgenome) >>syri1.sp1.sp2.chrnames;
  echo "###QRY-names###" >>syri1.sp1.sp2.chrnames;
  paste <(grep '^>' querygenlisted1.fasta) <(grep '^>' qrygenome) >>syri1.sp1.sp2.chrnames;       # CREATE pair aliases
  printf "\tExecuting minimap for 1st syri iteration - Second\n";
  minimap2 -t $(eval nproc) -k 21 -ax asm5 --score-N 2 -N 2 -p 0.9 -r 1000 -g 5000 --cs=long -Lc refgenome qrygenome | awk -F'\t' '$1~"^@" || $5>=20' >align.sp1.sp2.sam;
  printf "\tConverting it to bam\n";        
  rm -f align.sp1.sp2.bam.tmp.* ;
  bash ~/RGP/scripts/sam2bam_withmatch.sh align.sp1.sp2.sam align.sp1.sp2.bam ; 
fi

printf "Running 1st iter of SYRI\n";
if [[ ! -f "log_syri1_rcnames" ]]; then
  rm -f align.sp1.sp2.sam ;     
  printf "\tExecuting SYRI - Initial\n";
  printf "#!/usr/bin/bash\n\n" >iter1.syri.run.bash ;
  echo "python3 $PATH_TO_SYRI -c align.sp1.sp2.bam -r refgenome -q qrygenome -k -F B --no-chrmatch --log WARN --novcf --nosr --nosv --nosnp --nc $(eval nproc) --prefix iter1.syri.sp1.sp2. 1>log_syri1 2>&1" >>iter1.syri.run.bash ;
  bash iter1.syri.run.bash; 
  samtools faidx qrygenome ;    
  grep -e 'has high fraction of inverted alignments' log_syri1 |
      sed -e 's/.* query genome (//' -e 's/)\. Ensure that.*//' >log_syri1_rcnames;    # GET probable Rev-Compl
fi      

printf "RC if needed by large inverted seqs\n";
  count_rcnames=$(cat log_syri1_rcnames | wc -l) ;
  if [ $count_rcnames -gt 0 ]; then
    printf "" >qrygenome.listedby.sp1 ;
    printf "\tReordering the query genome\n";
    bash ~/RGP/scripts/rearrangements/convert_rc_4mlog-qrygen.sh log_syri1_rcnames qrygenome qrygenome.listedby.sp1 ; 
  else
    ln -sf qrygenome qrygenome.listedby.sp1 ;
  fi

printf "Creating qry seq as needed after 1st iter\n";  
if [[ ! -f "qrygenome.listedby.sp1.edit.fai" ]]; then
  samtools faidx qrygenome.listedby.sp1 ;
  cut -f1 qrygenome.listedby.sp1.fai >qrygenome.listedby.sp1.names ;               # GET names only for RC
  cut -f1 qrygenome.fai |
      grep -w -v -f qrygenome.listedby.sp1.names |
      xargs samtools faidx qrygenome >>qrygenome.listedby.sp1 ;                      # ADD rest of the sequences
  samtools faidx qrygenome.listedby.sp1 ;
  samtools faidx qrygenome ; 
  awk -F'\t' '{print $1}' qrygenome.fai >qrygenome.fai.names; 
  cat qrygenome.fai.names | xargs samtools faidx qrygenome.listedby.sp1 >qrygenome.listedby.sp1.edit ;
  samtools faidx qrygenome.listedby.sp1.edit ;
fi

  printf "\tCreating intermediate files before SYRI\n";
  bash ~/RGP/scripts/edit_chrheader.sh qrygenome.listedby.sp1.edit qrygenome_fin ;
  echo "###REF-names###" >syri2.sp1.sp2.chrnames;
  paste <(grep '^>' ${reference}) <(grep '^>' refgenome) >>syri2.sp1.sp2.chrnames;
  echo "###QRY-names###" >>syri2.sp1.sp2.chrnames;
  paste <(grep '^>' ${query}) <(grep '^>' qrygenome_fin) >>syri2.sp1.sp2.chrnames;  # CREATE PAIR aliases
  printf "\tExecuting minimap for SYRI - THIRD \n";
  rm -f align.sp1.sp2.bam ;
  minimap2 -t $(eval nproc) -ax asm10 --score-N 2 -g 1000 --cs=long -Lc refgenome qrygenome_fin | awk -F'\t' '$1~"^@" || $5>=20' >align.sp1.sp2.sam;
  printf "\tBAMing for SYRI - 2nd iteration\n";        
  bash ~/RGP/scripts/sam2bam_withmatch.sh align.sp1.sp2.sam bam ; 
  rm -f align.sp1.sp2.sam ;
  printf "\tExecute SYRI - Part II \n" ;

  module load minimap2 samtools bedtools ;
  printf "\tExecuting SYRI - FINALLY\n" ;
  python3 $PATH_TO_SYRI -c bam -r refgenome -q qrygenome_fin -b 5 -k -F B --no-chrmatch --nosnp --lf Final.syri.LOG --nc $(eval nproc) --cigar --prefix Final.sp1.sp2. 1>LOG_syri2 2>&1;
  bash ~/RGP/scripts/rearrangements/execute_synt-replacer.sh Final.sp1.sp2.syri.out syri1.sp1.sp2.chrnames Final.sp1.sp2.syri.synt.list ;
  bash ~/RGP/scripts/rearrangements/execute_syriout-namereplacer.sh Final.sp1.sp2.syri.out syri1.sp1.sp2.chrnames Final.sp1.sp2.syri.out.TSV ;

