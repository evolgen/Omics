#!/usr/bin/bash

mkdir -p new_Elements new_Scores            # put fragments here

file=$@

#for file in Chunks/*.ss
#do
	beds=$(echo $file | sed -e 's/new_Chunks/new_Elements/' -e 's/\.ss/.bed/')
	scores=$(echo $file | sed -e 's/new_Chunks/new_Scores/' -e 's/\.ss//')
	phastCons --gc 0.4 --most-conserved ${beds}.1 --score $file Lacerta_ave.cons.mod_1,Lacerta_ave.noncons.mod_1 >${scores}.1
	phastCons --gc 0.4 --most-conserved ${beds}.2 --score $file Lacerta_ave.cons.mod_2,Lacerta_ave.noncons.mod_2 >${scores}.2
	phastCons --gc 0.4 --most-conserved ${beds}.3 --score $file Lacerta_ave.cons.mod_3,Lacerta_ave.noncons.mod_3 >${scores}.3
	phastCons --gc 0.4 --most-conserved ${beds}.4 --score $file Lacerta_ave.cons.mod_4,Lacerta_ave.noncons.mod_4 >${scores}.4
	phastCons --gc 0.4 --most-conserved ${beds}.5 --target-coverage 0.1 --expected-length 8 --score $file Lacerta_ave.cons.mod_5,Lacerta_ave.noncons.mod_5 >${scores}.5
#	phastCons --gc 0.4 --most-conserved ${beds}.6 --target-coverage 0.1 --expected-length 8 --score $file Lacerta_ave.cons.mod_6,Lacerta_ave.noncons.mod_6 >${scores}.6
	echo "$file	DONE"
	
#done



