
SV<-matrix(c(2363,3004,121,787),nrow=2, byrow="T")   #0
fisher.test(SV,alternative="two.sided")
fisher.test(SV,alternative="greater")
fisher.test(SV,alternative="less")

complex      positive        0       503     0               6389
complex      positive        0       5367     0               1233
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0 Inf
# sample estimates:
#   odds ratio 
# 0 

DEL  positive        140     363     2329            3638
# p-value = 2.68e-07
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.7164918
# sample estimates:
#   odds ratio 
# 0.6024871 

DEL  neutral 2363    3004    121             787
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   4.316097      Inf
# sample estimates:
#   odds ratio 
# 5.114828 

DUP  positive        20      483     211             6138
# p-value = 0.8202
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.000000 1.804566
# sample estimates:
#   odds ratio 
# 1.204521 

DUP  neutral 215     5152    17              1189
# p-value = 9.665e-07
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   1.900228      Inf
# sample estimates:
#   odds ratio 
# 2.918447 

IDP  positive        0       503     2               6387
IDP  neutral 2       5365    0               1233
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.04313348        Inf
# sample estimates:
#   odds ratio 
# Inf 

INS  positive        84      419     1481            4657
# p-value = 6.07e-05
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.7755632
# sample estimates:
#   odds ratio 
# 0.6304198 

INS  neutral 1495    3872    74              969
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   4.103857      Inf
# sample estimates:
#   odds ratio 
# 5.05489 

INVDUP       positive        0       503     0               6389
INVDUP       neutral 0       5367    0               1233
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0 Inf
# sample estimates:
#   odds ratio 
# 0 

INV  positive        10      493     135             6216
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4350419 1.7864273
# sample estimates:
#   odds ratio 
# 0.9339709 

INV  neutral 139     5228    6               1192
# p-value = 2.91e-07
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   2.612919      Inf
# sample estimates:
#   odds ratio 
# 5.281309 

TRA  positive        0       503     2               6386
TRA  neutral 2       5365    0               1232
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.04309852        Inf
# sample estimates:
#   odds ratio 
# Inf 


SV<-matrix(c(2,5365,0,1232),nrow=2, byrow="T")   #0
fisher.test(SV,alternative="two.sided")
fisher.test(SV,alternative="greater")
fisher.test(SV,alternative="less")
