#!/usr/bin/bash

mkdir -p Elements Scores            # put fragments here

file=$@

#for file in Chunks/*.ss
#do
	beds=$(echo $file | sed -e 's/Chunks/Elements/' -e 's/\.ss/.bed/')
	scores=$(echo $file | sed -e 's/Chunks/Scores/' -e 's/\.ss//')
	phastCons --gc 0.4 --most-conserved ${beds}.1 --score $file ave.cons.mod_1,ave.noncons.mod_1 >${scores}.1
	phastCons --gc 0.4 --most-conserved ${beds}.2 --score $file ave.cons.mod_2,ave.noncons.mod_2 >${scores}.2
	phastCons --gc 0.4 --most-conserved ${beds}.3 --score $file ave.cons.mod_3,ave.noncons.mod_3 >${scores}.3
	phastCons --gc 0.4 --most-conserved ${beds}.4 --score $file ave.cons.mod_4,ave.noncons.mod_4 >${scores}.4
	phastCons --gc 0.4 --most-conserved ${beds}.5 --score $file ave.cons.mod_5,ave.noncons.mod_5 >${scores}.5
	echo "$file	DONE"
	
#done



