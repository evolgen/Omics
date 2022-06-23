
# 33824 ../annotation/Lacerta_viridis_annotation_FINAL.bed
# 7917  Total DE-genes


LVB_Pos_SVs<-matrix(c(7,129,140,1801),nrow=2, byrow="T")
#LVB_Pos_SVs<-matrix(c(41,2934,30,2715),nrow=2, byrow="T")
fisher.test(LVB_Pos_SVs,alternative="less")
### Positive results


1) W_CpG_SVs<-matrix(c(44,8627,30,8696),nrow=2, byrow="T")
#LVB_Pos_SVs<-matrix(c(41,2934,30,2715),nrow=2, byrow="T")
# fisher.test(W_CpG_SVs,alternative="greater")
# p-value = 0.06131
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.976786      Inf
# sample estimates:
#   odds ratio 
# 1.478357 

1) LVB_Pos_SVs<-matrix(c(41,2934,30,2715),nrow=2, byrow="T")
#LVB_Pos_SVs<-matrix(c(41,2934,30,2715),nrow=2, byrow="T")
# p-value = 0.3417
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7681349 2.1038755
# sample estimates:
#   odds ratio 
# 1.264602 

1) All_SVS <- 198,9791,191,4087
# p-value = 8.989e-16
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.5155054
# sample estimates:
#   odds ratio 
# 0.4327467 


2) DEL_CDS <- matrix(c(2,123,12,78),nrow=2, byrow="T")
p-value = 0.0006962
alternative hypothesis: true odds ratio is less than 1
95 percent confidence interval:
0.000000 0.415164
sample estimates:
odds ratio 
0.1067202 


#fisher.test(All_SVs,alternative="two.sided")
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5761904 0.6399852
# sample estimates:
#   odds ratio 
# 0.6073011 



# Sanity test
# All_SVs<-matrix(c(2797,12264,5120,13643),nrow=2, byrow="T")   #0
# fisher.test(All_SVs,alternative="greater")

#             P-Sel   !P-Sel
# SV
# Other-SV


# [rohit@bloodymary SVs]$ bedtools intersect -u -a <(bedtools intersect -v -a SV_Deletion.bed -b <(cat SV_Duplication.bed SV_Insertion.bed SV_Inversion.bed SV_Unknown.bed SV_Transposition.bed)) -b <(awk -F'\t' '$4<=1' Selection_coordinates_CDS.bed_nocompl_DelTraUnk_filt | bedtools intersect -u -a /scr/bloodymary/rohit/Lacerta_viridis/annotation/Lacerta_viridis_cdsbasedannotation.CDS.bed -b -) | wc -l
# 123
# [rohit@bloodymary SVs]$ bedtools intersect -u -a <(bedtools intersect -v -a SV_Deletion.bed -b <(cat SV_Duplication.bed SV_Insertion.bed SV_Inversion.bed SV_Unknown.bed SV_Transposition.bed)) -b <(awk -F'\t' '$4>1' Selection_coordinates_CDS.bed_nocompl_DelTraUnk_filt | bedtools intersect -u -a /scr/bloodymary/rohit/Lacerta_viridis/annotation/Lacerta_viridis_cdsbasedannotation.CDS.bed -b -) | wc -l
# 2
# [rohit@bloodymary SVs]$ bedtools intersect -u -a <(bedtools intersect -v -b SV_Deletion.bed -a <(cat SV_Duplication.bed SV_Insertion.bed SV_Inversion.bed SV_Unknown.bed SV_Transposition.bed)) -b <(awk -F'\t' '$4<=1' Selection_coordinates_CDS.bed_nocompl_DelTraUnk_filt | bedtools intersect -u -a /scr/bloodymary/rohit/Lacerta_viridis/annotation/Lacerta_viridis_cdsbasedannotation.CDS.bed -b -) | wc -l
# 78
# [rohit@bloodymary SVs]$ bedtools intersect -u -a <(bedtools intersect -v -b SV_Deletion.bed -a <(cat SV_Duplication.bed SV_Insertion.bed SV_Inversion.bed SV_Unknown.bed SV_Transposition.bed)) -b <(awk -F'\t' '$4>1' Selection_coordinates_CDS.bed_nocompl_DelTraUnk_filt) | bedtools intersect -u -a /scr/bloodymary/rohit/Lacerta_viridis/annotation/Lacerta_viridis_cdsbasedannotation.CDS.bed -b - | wc -l
# 12
