
INV<-matrix(c(602,3,32753,277),nrow=2, byrow="T")   #0
fisher.test(INV,alternative="two.sided")
fisher.test(INV,alternative="greater")
fisher.test(INV,alternative="less")

del_H<-matrix(c(2133,9994,3150,18547),nrow=2, byrow="T")   #12127
fisher.test(del_H,alternative="two.sided")
fisher.test(del_H,alternative="less")
fisher.test(del_H,alternative="greater")
# p-value = 7.059e-14
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   1.194194      Inf
# sample estimates:
#   odds ratio 
# 1.256619 
p_del_H <- fisher.test(del_H,alternative="greater")$p.value
p.adjust(p_del_H, method = p.adjust.methods, n = length(p_del_H))
#[1] 7.058765e-14


dup_H<-matrix(c(141,735,5142,27806),nrow=2, byrow="T")   #876
fisher.test(dup_H,alternative="greater")
# p-value = 0.361
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.8845396       Inf
# sample estimates:
#   odds ratio 
# 1.037385 

inv_H<-matrix(c(163,769,5120,27772),nrow=2, byrow="T")   #932
fisher.test(inv_H,alternative="greater")
# p-value = 0.06226
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.9902382       Inf
# sample estimates:
#   odds ratio 
# 1.149696 
p_inv_H <- fisher.test(inv_H,alternative="greater")$p.value
p.adjust(p_inv_H, method = p.adjust.methods, n = length(p_inv_H))
#[1] 0.06225634


del_L<-matrix(c(2122,10005,4339,17496),nrow=2, byrow="T")   #12127
fisher.test(del_L,alternative="greater")
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.8145963       Inf
# sample estimates:
#   odds ratio 
# 0.8552311 

dup_L<-matrix(c(138,738,6323,26625),nrow=2, byrow="T")   #876
fisher.test(dup_L,alternative="greater")
# p-value = 0.996
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.6706192       Inf
# sample estimates:
#   odds ratio 
# 0.7873961 

inv_L<-matrix(c(161,771,6300,26592),nrow=2, byrow="T")   #932
fisher.test(inv_L,alternative="greater")
# p-value = 0.9321
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.7588386       Inf
# sample estimates:
#   odds ratio 
# 0.8814217 


del_DE<-matrix(c(4255,7872,7489,14208),nrow=2, byrow="T")   #12127
fisher.test(del_DE,alternative="greater")
# p-value = 0.1479
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.9858274       Inf
# sample estimates:
#   odds ratio 
# 1.025485 

dup_DE<-matrix(c(279,597,11465,21483),nrow=2, byrow="T")   #876
fisher.test(dup_DE,alternative="greater")
# p-value = 0.9681
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.7735294       Inf
# sample estimates:
#   odds ratio 
# 0.8756789 

inv_DE<-matrix(c(324,608,11420,21472),nrow=2, byrow="T")   #932
fisher.test(inv_DE,alternative="greater")
# p-value = 0.5014
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.8907947       Inf
# sample estimates:
#   odds ratio 
# 1.001954 

idp_DE<-matrix(c(1,1,11743,22079),nrow=2, byrow="T")   #2
fisher.test(idp_DE,alternative="greater")
# p-value = 0.5739
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   0.04883818        Inf
# sample estimates:
#   odds ratio 
# 1.880216 

total_SVs <- c("132659")
total_genes_withSVs <- c("16361")
total_diff_genes <- c("280")

total_genes_withSVs_new <- c("13937")
total_diff_genes_new <- c("11744")
total_genes_new <- c("33824")


# view the matrix

# fischer test with different options
fisher.test(m,alternative="two.sided")
fisher.test(m,alternative="greater")
fisher.test(m,alternative="less")