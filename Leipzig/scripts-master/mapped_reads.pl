cut --fields=1 ./tophat_out/accepted_hits.sam | sort --unique > ./tophat_out/accepted_readIds.txt
paste - - - - < ./reads_1.fq | fgrep --invert-match --file ./tophat_out/accepted_readsIds.txt | tr "\t" "\n" > ./readsFilt_1.fq
