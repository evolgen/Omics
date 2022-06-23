#!/usr/bin/bash

#mkdir -p Chunks Trees            # put fragments here

file=$@

#for file in Chunks/*.ss
#do
	root2=$(echo $file | sed -e 's/new_Chunks/new_Trees/' -e 's/\.ss//')
	phastCons --gc 0.4 --estimate-trees $root2 $file Lacerta_init.mod --no-post-probs 2>>Lacerta_check_phastcons_2
	echo "$root	DONE"
	
#done



