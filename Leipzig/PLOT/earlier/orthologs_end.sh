#!/usr/bin/bash


#Remove the redundant unitigs using the VMatch tool

#1) First, go for blat on the two genomes
pblat -threads=8 -noHead target_unitigs.fa reference_unitigs.fa blat1.psl
pblat -out=blast8 -threads=8  reference_unitigs.fa target_unitigs.fa blat2t.blast8

#2) Filter by evalue of blast-like format
awk '$11<1e-15' blat2t.blast8 > blat2.blast8
rm blat2t.blast8

#3) Extract orthologous hits by RBH ####No-hits => missing data or unique regions
awk -F'	' 'NR==FNR{c[$10$14]++;next};c[$2$1] > 0' blat2.blast8 blat1.psl > 1st_hits
rm blat2.blast8 blat1.psl
awk '{print $10}' 1st_hits | sort -u > ref_names
awk '{print $14}' 1st_hits | sort -u > tar_names

#4) Use samtools to extract these sequences
samtools faidx reference_unitigs.fa
xargs samtools faidx reference_unitigs.fa < ref_names > reference_unitigs_filt.fa
samtools faidx target_unitigs.fa
xargs samtools faidx target_unitigs.fa < tar_names > target_unitigs_filt.fa

#5)  Use LastZ to align sequences properly as it uses a better scoring matrix, convert maf to psl format
#lastz target_unitigs_filt.fa[multiple] --step=18 reference_unitigs_filt.fa --format=maf >lastz_hits.maf
lastz target_unitigs_filt.fa[multiple] --step=9 --noytrim --chain --identity=0..100 --ambiguous=iupac --hspthresh=6000 ‑‑traceback=2000.0M --inner=2500 --masking=254 reference_unitigs_filt.fa --format=maf >lastz_hits.maf
maf-convert psl lastz_hits.maf >lastz_hits.psl

#6) Extract bed cordinates from the lastz psl outfile
perl blat_insbl_ext_bed_all.pl lastz_hits.psl

#7) Extract the ones longer the length-threshold
awk '$7>=20' lastz_hits_query.bed > query.bed
awk '$7>=20' lastz_hits_target.bed > target.bed
awk '$7<17' lastz_hits_query.bed > small_query.bed
awk '$7<17' lastz_hits_target.bed > small_target.bad



#8) Seqtk to extract the sequences depending on the customized bed-format
seqtk subseq reference_unitigs_filt.fa query.bed > query_bed.fasta
seqtk subseq target_unitigs_filt.fa target.bed > target_bed.fasta 

#9) Blast search for highly similar sequences
makeblastdb -in reference_unitigs.fa -dbtype 'nucl' -out ref_db
makeblastdb -in target_unitigs.fa dbtype 'nucl' -out tar_db
blastn -task megablast -query query_bed.fasta -db ref_db -out query_ins_self_hits -perc_identity 90 -word_size 16 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"
awk -F'	' 'NR==FNR{c[$1$2]++;next};c[$10$14] > 0' lastz_hits.psl query_ins_self_hits > query_ins_self_hits_ortholog                               #Inversion
blastn -task megablast -query target_bed.fasta -db tar_db -out target_ins_self_hits -perc_identity 90 -word_size 16 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"
awk -F'	' 'NR==FNR{c[$1$2]++;next};c[$10$14] > 0' lastz_hits.psl target_ins_self_hits > target_ins_self_hits_ortholog                                #Inversion
blastn -task megablast -query query_bed.fasta -db tar_db -out query_ins_target_hits -perc_identity 90 -word_size 16 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"
#awk -F'	' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' query_ins_target_hits lastz_hits.psl > 1st_hits
blastn -task megablast -query target_bed.fasta -db ref_db -out target_ins_query_hits -perc_identity 90 -word_size 16 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"
#awk -F'	' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' target_ins_query_hits lastz_hits.psl > 1st_hits



######### EDIT





[rohit@K53 real]$ sed 's/:.*//p' que_ins_vcrbd_hits > qh1; cut -f2 que_ins_vcrbd_hits > qh2; paste qh1 qh2 > qh; rm qh1 qh2
[rohit@K53 real]$ awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' qh lastz_accurate_orthologs_ranked_indels.psl | head

[rohit@K53 real]$ sed 's/:.*//p' tar_ins_ccrbd_hits > th1; cut -f2 tar_ins_ccrbd_hits > th2; paste th1 th2 >th; rm th1 th2
[rohit@K53 real]$ awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$2$1] > 0' th lastz_accurate_orthologs_ranked_indels.psl



#####For filter file
Extract query and target beds
time perl blat_insbl_ext_bed_all.pl lastz_accurate_orthologs_ranked_indels_filter.psl
Extract the ones larger than 17 bases only for query and target, with direction for query
[rohit@K53 filter]$ awk 'BEGIN {OFS="\t"} {if( $7>17 && $6=="+")print $0;}' lastz_accurate_orthologs_ranked_indels_filter_query_all.bed >lastz_accurate_orthologs_ranked_indels_filter_query_all_big_F.bed; awk 'BEGIN {OFS="\t"} {if( $7>17 && $6=="-")print $0;}' lastz_accurate_orthologs_ranked_indels_filter_query_all.bed >lastz_accurate_orthologs_ranked_indels_filter_query_all_big_R.bed; awk 'BEGIN {OFS="\t"} {if( $7>17 && $6=="+")print $0;}' lastz_accurate_orthologs_ranked_indels_filter_target_all.bed >lastz_accurate_orthologs_ranked_indels_filter_target_all_big_F.bed

Extract all the sequences
[rohit@K53 real]$ /scr/bloodymary/rohit/seqtk-master/seqtk subseq clint_cr_seq_clean.fa lastz_accurate_orthologs_ranked_indels_filter_query_all_big_F.bed > lastz_accurate_orthologs_ranked_indels_filter_query_all_big_F.fasta; /scr/bloodymary/rohit/seqtk-master/seqtk subseq clint_cr_seq_clean.fa lastz_accurate_orthologs_ranked_indels_filter_query_all_big_R.bed > lastz_accurate_orthologs_ranked_indels_filter_query_all_big_R.fasta; /scr/bloodymary/rohit/seqtk-master/seqtk subseq vail_cr_seq_clean.fa lastz_accurate_orthologs_ranked_indels_filter_target_all_big_F.bed > lastz_accurate_orthologs_ranked_indels_filter_target_all_big_F.fasta 

[rohit@K53 filter]$  time blastn -task megablast -query lastz_accurate_orthologs_ranked_indels_filter_query_all_big_F.fasta -db ../v_cr_db -out que_F_ins_vcrbd_hits -perc_identity 98 -word_size 16 -evalue 1e-50 -num_threads 6 -max_target_seqs 3 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"; time blastn -task megablast -query lastz_accurate_orthologs_ranked_indels_filter_query_all_big_R.fasta -db ../v_cr_db -out que_R_ins_vcrbd_hits -perc_identity 98 -word_size 16 -evalue 1e-50 -num_threads 6 -max_target_seqs 3 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"; time blastn -task megablast -query lastz_accurate_orthologs_ranked_indels_filter_target_all_big_F.fasta -db ../c_cr_db -out tar_F_ins_ccrbd_hits -perc_identity 98 -evalue 1e-50 -num_threads 6 -max_target_seqs 3 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"

[rohit@K53 filter]$ sed 's/:.*//p' que_F_ins_vcrbd_hits > qh1; cut -f2 que_F_ins_vcrbd_hits > qh2; paste qh1 qh2 > qhF; rm qh1 qh2; awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' qhF lastz_accurate_orthologs_ranked_indels_filter.psl | wc -l
0
[rohit@K53 filter]$ sed 's/:.*//p' que_R_ins_vcrbd_hits > qh1; cut -f2 que_R_ins_vcrbd_hits > qh2; paste qh1 qh2 > qhR; rm qh1 qh2; awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' qhR lastz_accurate_orthologs_ranked_indels_filter.psl | wc -l
0
[rohit@K53 filter]$ sed 's/:.*//p' tar_F_ins_ccrbd_hits > qh1; cut -f2 tar_F_ins_ccrbd_hits > qh2; paste qh1 qh2 > thF; rm qh1 qh2; awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' thF lastz_accurate_orthologs_ranked_indels_filter.psl | wc -l
0

Check how many ins_gaps got hits with only CR sequences
[rohit@K53 filter]$ sort -k1,1 que_F_ins_vcrbd_hits | uniq | wc -l
11
[rohit@K53 filter]$ sort -k1,1 que_R_ins_vcrbd_hits | uniq | wc -l
9
[rohit@K53 filter]$ sort -k1,1 tar_F_ins_ccrbd_hits | uniq | wc -l
38


Compare with all contigs
[rohit@K53 filter]$ makeblastdb -dbtype nucl -in ../../../clint_79_k29_ok.fa -out clint_db; makeblastdb -dbtype nucl -in ../../../vailant_58_k27_ok.fa -out vailant_db
time blastn -task megablast -query lastz_accurate_orthologs_ranked_indels_filter_query_all_big_F.fasta -db ../vailant_db -out all_que_F_ins_vcrbd_hits -word_size 16 -evalue 1e-50 -num_threads 6 -max_target_seqs 3 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"; time blastn -task megablast -query lastz_accurate_orthologs_ranked_indels_filter_query_all_big_R.fasta -db ../vailant_db -out all_que_R_ins_vcrbd_hits -word_size 16 -evalue 1e-50 -num_threads 6 -max_target_seqs 3 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"; time blastn -task megablast -query lastz_accurate_orthologs_ranked_indels_filter_target_all_big_F.fasta -db ../clint_db -out all_tar_F_ins_ccrbd_hits -evalue 1e-50 -num_threads 6 -max_target_seqs 3 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand"

Check how many ins_gaps got hits with all assembly
[rohit@K53 filter]$ sort -k1,1 all_que_F_ins_vcrbd_hits | uniq | wc -l;  sort -k1,1 all_que_R_ins_vcrbd_hits | uniq | wc -l; sort -k1,1 all_tar_F_ins_ccrbd_hits | uniq | wc -l
23
15
68

[rohit@K53 filter]$ sed 's/:.*//p' all_que_F_ins_vcrbd_hits > qh1; cut -f2 all_que_F_ins_vcrbd_hits > qh2; paste qh1 qh2 > qhF; rm qh1 qh2; awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' qhF ../lastz_accurate_orthologs_ranked_indels_filter.psl | wc -l
0
[rohit@K53 filter]$ sed 's/:.*//p' all_que_R_ins_vcrbd_hits > qh1; cut -f2 all_que_R_ins_vcrbd_hits > qh2; paste qh1 qh2 > qhR; rm qh1 qh2; awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' qhR lastz_accurate_orthologs_ranked_indels_filter.psl | wc -l
0
[rohit@K53 filter]$ sed 's/:.*//p' all_tar_F_ins_ccrbd_hits > th1; cut -f2 all_tar_F_ins_ccrbd_hits > th2; paste th1 th2 > thF; rm th1 th2; awk -F'|' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' thF lastz_accurate_orthologs_ranked_indels_filter.psl | wc -l
0



e-value for the lastz accurate ortholog filter hits
[rohit@K53 filter]$ awk -F'     ' 'NR==FNR{c[$1$2]++;next};c[$10$14] > 0' lastz_accurate_orthologs_ranked_indels_filter.psl ../vaildb_clint.blast8  > eval_lastz_hits
#[rohit@K53 filter]$ awk -F'     ' 'NR==FNR{c[$10$14]++;next};c[$1$2] > 0' ../vaildb_clint.blast8 lastz_accurate_orthologs_ranked_indels_filter.psl > eval_lastz_hits_2                    
[rohit@K53 filter]$ sort -k1,1 -k2,2 eval_lastz_hits | awk '!x[$1,$2]++' | wc -l
1966981
[rohit@K53 filter]$ wc -l eval_lastz_hits
2410829 eval_lastz_hits
#[rohit@K53 filter]$ sort -k10,10 -k14,14 eval_lastz_hits_2 | awk '!x[$10,$14]++' | wc -l
#42432
#[rohit@K53 filter]$ wc -l eval_lastz_hits_2
#42639 eval_lastz_hits_2


################ Tests
Jim kent's advice
[rohit@K53 real]$ pslPretty accurate_orthologs_ranked_indels vail_cr_seq_clean.fa clint_cr_seq_clean.fa accurate_orthologs_ranked_indels_prettypsl.axt
Run threshold on bloodymary (32-bit only)
[rohit@bloodymary real]$ /opt/kentUtils-master/bin/subsetAxt accurate_orthologs_ranked_indels_prettypsl.axt accurate_orthologs_ranked_indels_prettypsl_threshold.axt /opt/kentUtils-master/src/hg/mouseStuff/subsetAxt/90.mat 3000
This extracts only rhe alifnment bloacks which cross the threshold. We get the results as bed files with line numbers, need to edit them - The problem is that we miss the high similar ones with wrong gaps as they do not cross the threshold


Extract the seqeunces for the insertion blocks so that they can be searched later
[rohit@K53 real]$ perl blat_insbl_ext.pl accurate_orthologs_ranked

Take the seqeunce files for which we found blat hits in case of indels and use them for extracting the sub-sequences
[rohit@K53 real]$ cp other_indel/*cr_seq_clean* .





Combine both small seqeunces files, and then go for clustering
[rohit@K53 real]$ bedtools getfasta -s -fi clint_cr_seq_clean.fa -name -bed accurate_orthologs_ranked_query_small.bed -s -fo accurate_orthologs_ranked_query_insert_extract_small.fa; bedtools getfasta -s -fi vail_cr_seq_clean.fa -name -bed accurate_orthologs_ranked_target_small.bed -s -fo accurate_orthologs_ranked_target_insert_extract_small.fa
Make a reference file for cluster representation
[rohit@K53 real]$ cat accurate_orthologs_ranked_query_insert_extract_small.fa accurate_orthologs_ranked_target_insert_extract_small.fa > accurate_orthologs_ranked_insert_extract_small_qt.fa
Cluster the small seqeunces so that they have maximum length covered, (in a cluster smaller seqeunces are completely covered), go for the best cluster than the first hit
[rohit@K53 real]$ time cd-hit-est -i accurate_orthologs_ranked_insert_extract_small_qt.fa -o qt_small_clust -n 4 -l 5 -M 20000 -d 0 -s 0.9 -aL 0.9 -aS 1 -g 1
Make two clusters, one with high coverage for smaler seqeunce (_high), one with no such criterion (_low)
[rohit@K53 real]$ time cd-hit-est -i accurate_orthologs_ranked_insert_extract_small_qt.fa -o qt_small_clust_high -n 4 -l 5 -M 20000 -d 0 -s 0.9 -aS 0.9 -g 1; time cd-hit-est -i accurate_orthologs_ranked_insert_extract_small_qt.fa -o qt_small_clust_low -n 4 -l 5 -M 20000 -d 0 -aS 0.9 -g 1
Cluster-check
[rohit@K53 real]$ grep -B 1 'Cluster ' qt_small_clust.clstr | grep -v 'Cluster' | awk '$1>0' | tail


[rohit@K53 real]$ bedtools getfasta -s -fi clint_cr_seq_clean.fa -name -bed accurate_orthologs_ranked_query.bed -s -fo accurate_orthologs_ranked_query_insert_extract.fa; bedtools getfasta -s -fi vail_cr_seq_clean.fa -name -bed accurate_orthologs_ranked_target.bed -s -fo accurate_orthologs_ranked_target_insert_extract.fa
[rohit@K53 real]$ time blastn -task megablast -query accurate_orthologs_ranked_query_insert_extract.fa -db v_cr_db -out que_ins_vcrbd_hits -perc_identity 95 -word_size 16 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt 6; time blastn -task megablast -query accurate_orthologs_ranked_target_insert_extract.fa -db c_cr_db -out tar_ins_ccrbd_hits -perc_identity 95 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt 6




[rohit@K53 real]$ makeblastdb -dbtype nucl -in accurate_orthologs_ranked_target_insert_extract.fa -out target_ins_db; makeblastdb -dbtype nucl -in accurate_orthologs_ranked_query_insert_extract.fa -out query_ins_db

[rohit@K53 real]$ time blastn -task megablast -query accurate_orthologs_ranked_query_insert_extract.fa -db query_ins_db -out tar_ins_qryinsdb_hits -perc_identity 95 -word_size 16 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt 6; time blastn -task blastn-short -query accurate_orthologs_ranked_target_insert_extract.fa -db target_ins_db -out que_ins_tarinsdb_hits -perc_identity 95 -evalue 1e-20 -num_threads 6 -max_target_seqs 5 -outfmt 6
############


##################	Finding Definite Indels	######################	4
Make seperate folders for Indels/Translocations and Inversions
Group into 4 categories based on sizes of insertion block
			1) <=5			2) 6-16			3) 17-30			4) >30
Make folders for each size - [rohit@k50 real]$ mkdir 5bp 16bp 30bp 100bp
1) <=5
Can be due to mismatches or polymorphism
But are considered indels, as we check only the insertion blocks

Find if multiple or single crs and others
[rohit@k50 real]$ perl blat_threshold_indel.pl accurate_orthologs_ranked
Move the indels, other SVs, snps, perfect orthologs with no gaps  to other_indels folder
      0 accurate_orthologs_ranked_multiple_cr_inone
  10692 accurate_orthologs_ranked_multiple_exact_cr
   3226 accurate_orthologs_ranked_multiple_uneven_cr
      0 accurate_orthologs_ranked_one_idt
  11056 accurate_orthologs_ranked_remaining
  13291 accurate_orthologs_ranked_single_cr_both
  38265 total


Go to the folder of other indels, deal with each indel seperately


################	4


Use Reapr tool for finding the contigs that are not false overlaps
[rohit@k50 chimp]$ sed -e 's/.*REAPR stats.*vailant/vailant/g' -e 's/. Assuming.*//g' reapr_warning | sed "s/'\.//g" | sort > vail_bad_seqnames


Remove the PSL results of these bad contigs, use the rest for the analysis
[rohit@k50 chimp]$ #grep -v -f /scr/k48san/rohit/chimp/vail_bad_seqnames lastz_accurate_orthologs_ranked_indels.psl > lastz_accurate_orthologs_ranked_indels_filter.psl





############
1
#pslPretty
#rohit@k52 /scr/bloodymary/rohit/chimp/CRs/real $ pslPretty accurate_orthologs_ranked ../../vailant_58_k27_R_soapcontigs.fa ../../clint_79_k29_soapcontigs.fa accurate_orthologs_ranked_pretty.out
#pslReps
#[rohit@K53 real]$ pslReps accurate_orthologs_ranked accurate_orthologs_ranked_reps.psl accurate_orthologs_ranked_reps.psr
1

