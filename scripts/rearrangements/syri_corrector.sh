#!/usr/bin/bash

set -e 

\ls -d /global/scratch2/rohitkolora/Rockfish/Genomes/FREEZE/*/*/Pairwise/SYRI/*_*/Final.*.*.syri.out |
    while read alnout; do
        dir_name=$(dirname "$alnout") ;
        reference=$(echo $alnout | sed -e 's/\/Pairwise\/SYRI\/.*/\/long_seq.fa/') ;
        chrsyri2=$(echo $alnout | sed -e 's/\/Final/\/syri1./' -e 's/.syri.out$/.chrnames2/') ; 
        syntblocks=$(echo $alnout | sed -e 's/\.syri\.out$/.syri.synt.list/') ;
        fintsvout=$(echo $alnout | sed -e 's/\.syri\.out$/.syri.out.TSV/') ;

        cd ${dir_name} ;
        module load minimap2 samtools bedtools ;
        echo "###REF-names###" >${chrsyri2} ;
        paste <(grep '^>' ${reference}) <(grep '^>' refgenome) >>${chrsyri2} ;
        echo "###QRY-names###" >>${chrsyri2} ;
        paste <(grep '^>' querygenlisted1.fasta) <(grep '^>' qrygenome) >>${chrsyri2} ;      # CREATE pair aliases
        
        bash ~/RGP/scripts/rearrangements/execute_synt-replacer.sh ${alnout} ${chrsyri2} ${syntblocks} ;
        bash ~/RGP/scripts/rearrangements/execute_syriout-namereplacer.sh ${alnout} ${chrsyri2} ${fintsvout} ;
        
        echo "${syntblocks}" ;

    done    
