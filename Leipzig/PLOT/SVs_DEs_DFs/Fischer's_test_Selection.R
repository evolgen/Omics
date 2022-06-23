
# 33824 ../annotation/Lacerta_viridis_annotation_FINAL.bed
# 7917  Total DE-genes

y1 = 37
y2 = 155
n1 = 37+1430
n2 = 155+2150

Paml_Rearr<-matrix(c(2,7,1675,1459),nrow=2, byrow="T")   
fisher.test(Paml_Rearr,alternative="two.sided")



## Positive results
A)kaks_Rearr<-matrix(c(8,87,189,3557),nrow=2, byrow="T")   
fisher.test(kaks_Rearr,alternative="greater")
# p-value = 0.1121
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.8260386       Inf
# sample estimates:
#   odds ratio 
# 1.730268 

B)PAML_Lac_SVs<-matrix(c(3,3639,5,3130),nrow=2, byrow="T")   
fisher.test(PAML_Lac_SVs,alternative="two.sided")
# p-value = 0.4839
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.08008529 2.65520207
# sample estimates:
#   odds ratio 
# 0.5161313 


C)kaks_SVs<-matrix(c(57,2744,140,2849),nrow=2, byrow="T")   
fisher.test(kaks_SVs,alternative="less")
# p-value = 1.249e-08
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.5545399
# sample estimates:
#   odds ratio 
# 0.4227571 



I)PAML_SVs<-matrix(c(65,8103,44,13634),nrow=2, byrow="T")
# p-value = 1.943e-06
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 1.771891      Inf
# sample estimates:
# odds ratio 
# 2.485546 

II)PAML+SVs_DE<-matrix(c(24,2318,41,5785),nrow=2, byrow="T")
# p-value = 0.09234
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.9177674       Inf
# sample estimates:
#   odds ratio 
# 1.460816 



1) Tran_CDS<-matrix(c(3,32,7,264),nrow=2, byrow="T")   
# p-value = 0.0938
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.7495505       Inf
# sample estimates:
#   odds ratio 
# 3.513613 


2) ALL_Introns <- matrix(c(194,45,879,4893),nrow=2, byrow="T")   
# p-value < 2.2e-16
#   alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   17.99283      Inf
# sample estimates:
#   odds ratio 
# 23.9697   


### Lesser  
1) All_SVs <- matrix(c(198,9791,191,4087),nrow=2, byrow="T")
# p-value = 8.989e-16
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.5155054
# sample estimates:
#   odds ratio 
# 0.4327467 


2) DEL_CDS <- matrix(c(2,123,12,78),nrow=2, byrow="T")
# p-value = 0.0006962
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
# 0.000000 0.415164
# sample estimates:
# odds ratio 
# 0.1067202 

3) ALL_Introns <- matrix(c(194,9745,879,4893),nrow=2, byrow="T")   
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.1269428
# sample estimates:
#   odds ratio 
# 0.1108333 



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
