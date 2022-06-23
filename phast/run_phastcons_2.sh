#!/usr/bin/bash

#mkdir -p Chunks Trees            # put fragments here

file=$@

#for file in Chunks/*.ss
#do
	root2=$(echo $file | sed -e 's/Chunks/Trees/' -e 's/\.ss//')
	phastCons --gc 0.4 --estimate-trees $root2 $file init.mod --no-post-probs 2>>check_phastcons_2
	echo "$root	DONE"
	
#done



