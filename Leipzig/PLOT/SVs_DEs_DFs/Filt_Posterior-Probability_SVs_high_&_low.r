
setwd("/homes/biertank/rohit/Downloads/scripts/SVs_DEs_DFs")

# Libraries you need
library(ggplot2)

# Loading Data
#grfTable <- read.table("/homes/biertank/rohit/Downloads/scripts/SVs_DEs_DFs/in_Table",sep = "\t",header = TRUE)

# DATA (here I just test for one Domain IPR000003)
# n1 = 15061 #nrow(grfTable[which(grfTable$Human.Specific ==1),]) # totall # of Human-spec no mater they have the occu or not
# y1 = 2796 # nrow(grfTable[which(grfTable$IPR000003!=0 & grfTable$Human.Specific==1 ),])  # Human-spec Occu only 
# n2 = 5121 + 13642 # nrow(grfTable[which(grfTable$Human.Specific ==0),]) # totall # of not human-spec no matter they occur or not 
# y2 = 5121  # nrow(grfTable[which(grfTable$IPR000003!=0 & grfTable$Human.Specific==0),])  # not human-spec Occu 


# DE  = 7917
# Genes = 33824
# SVs = 48331
# High_DE = 3705
# Low_DE  = 4214

psi_1 = rbeta(10000000, y1+1, (n1-y1)+1)  
psi_2 = rbeta(10000000, y2+1, (n2-y2)+1)
diff_psi1_psi2 = psi_1 - psi_2
diff_quantiles = quantile(diff_psi1_psi2,c(0.005,0.025,0.5,0.975,0.995)) #,na.rm=TRUE)
print(diff_quantiles,digits=2)
print(mean(psi_1 > psi_2))
print(mean(diff_psi1_psi2))

complex	high	1	2287	9		17053
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00082 -0.00062  0.00016  0.00188  0.00269 
# > print(mean(psi_1 > psi_2))
# [1] 0.6201676
# > print(mean(diff_psi1_psi2))
# [1] 0.000287717
complex	low	0	129	10		19211
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00065 -0.00038  0.00475  0.02742  0.03931 
# > print(mean(psi_1 > psi_2))
# [1] 0.9285877
# > print(mean(diff_psi1_psi2))
# [1] 0.007062456
DEL	high	881	1407	5477		11585
# 0.5%  2.5%   50% 97.5% 99.5% 
# 0.037 0.043 0.064 0.085 0.092 
# > print(mean(psi_1 > psi_2))
# [1] 1
# > print(mean(diff_psi1_psi2))
# [1] 0.0641206
DEL	low	53	76	6305		12916
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0233  0.0013  0.0837  0.1695  0.1965 
# > print(mean(psi_1 > psi_2))
# [1] 0.9768456
# > print(mean(diff_psi1_psi2))
# [1] 0.084172
DUP	high	68	2220	514		16548
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00920 -0.00716 -0.00017  0.00774  0.01042 
# > print(mean(psi_1 > psi_2))
# [1] 0.4822282
# > print(mean(diff_psi1_psi2))
# [1] -5.002121e-05
DUP	low	9	120	573		18648
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00082  0.00756  0.04433  0.09752  0.11756 
# > print(mean(psi_1 > psi_2))
# [1] 0.9939731
# > print(mean(diff_psi1_psi2))
# [1] 0.04647536
IDP	high	1	2287	5		17057
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00050 -0.00033  0.00039  0.00210  0.00291 
# > print(mean(psi_1 > psi_2))
# [1] 0.8032618
# > print(mean(diff_psi1_psi2))
# [1] 0.0005218934
IDP	low	0	129	6		19215
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00040 -0.00017  0.00495  0.02762  0.03961 
# > print(mean(psi_1 > psi_2))
# [1] 0.9539018
# > print(mean(diff_psi1_psi2))
# [1] 0.007271248
INS	high	507	1781	3030		14032
# 0.5%  2.5%   50% 97.5% 99.5% 
# 0.021 0.026 0.044 0.062 0.068 
# > print(mean(psi_1 > psi_2))
# [1] 0.9999997
# > print(mean(diff_psi1_psi2))
# [1] 0.04420931
INS	low	33	96	3504		15717
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0139  0.0059  0.0760  0.1554  0.1818 
# > print(mean(psi_1 > psi_2))
# [1] 0.9838193
# > print(mean(diff_psi1_psi2))
# [1] 0.07719776
INVDUP	high	0	2288	0		17062
INVDUP	low	0	129	0		19221
INV	high	47	2241	359		16703
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00768 -0.00604 -0.00026  0.00646  0.00876 
# > print(mean(psi_1 > psi_2))
# [1] 0.4679522
# > print(mean(diff_psi1_psi2))
# [1] -0.0001350496
INV	low	8	121	398		18823
# 0.5%   2.5%    50%  97.5%  99.5% 
# 0.0037 0.0113 0.0458 0.0969 0.1164 
# > print(mean(psi_1 > psi_2))
# [1] 0.9982499
# > print(mean(diff_psi1_psi2))
# [1] 0.04794479
TRA	high	1	2287	13		17049
# 0.5%     2.5%      50%    97.5%    99.5% 
# -1.1e-03 -9.0e-04 -6.6e-05  1.7e-03  2.5e-03 
# > print(mean(psi_1 > psi_2))
# [1] 0.4557828
# > print(mean(diff_psi1_psi2))
# [1] 5.299865e-05
TRA	low	0	129	14		19207
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00089 -0.00060  0.00454  0.02721  0.03917 
# > print(mean(psi_1 > psi_2))
# [1] 0.9038493
# > print(mean(diff_psi1_psi2))
# [1] 0.006855037


y1 = 0
n1 = y1+129
y2 = 14
n2 = y2+19207
psi_1 = rbeta(10000000, y1+1, (n1-y1)+1)  
psi_2 = rbeta(10000000, y2+1, (n2-y2)+1)
diff_psi1_psi2 = psi_1 - psi_2
diff_quantiles = quantile(diff_psi1_psi2,c(0.005,0.025,0.5,0.975,0.995)) #,na.rm=TRUE)
print(diff_quantiles,digits=2)
print(mean(psi_1 > psi_2))
print(mean(diff_psi1_psi2))

#

#### SVs and no-SVs with DE and no-DE
y1 = 1
y2 = 3
n1 = y1+2287
n2 = y2+2689
  
complex	high	1	2287	3		2689
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00344 -0.00261 -0.00059  0.00130  0.00213 
# 0.2418642
# -0.0006122043
complex	low	3	2689	1		2287

DEL	high	881	1407	1110		1582
# 0.5%     2.5%      50%    97.5%    99.5% 
# -6.3e-02 -5.4e-02 -2.7e-02  2.8e-05  8.6e-03 
# 0.0251225
# -0.02724228
DEL	low	1110	1582	881		1407

DUP	high	68	2220	110		2582
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.02460 -0.02134 -0.01107 -0.00079  0.00249 
# 0.0174376
# -0.01107209
DUP	low	110	2582	68		2220

IDP	high	1	2287	1		2691
# 0.5%     2.5%      50%    97.5%    99.5% 
# -2.2e-03 -1.5e-03  9.9e-05  1.9e-03  2.7e-03 
# 0.5607045
# 0.0001310963
IDP	low	1	2691	1		2287

INS	high	507	1781	619		2073
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0388 -0.0315 -0.0083  0.0150  0.0224 
# 0.2418649
# -0.008307205
INS	low	619	2073	507		1781

INVDUP	high	0	2288	0		2692
INVDUP	low	0	2692	0		2288

INV	high	47	2241	66		2626
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0149 -0.0122 -0.0039  0.0044  0.0071 
# 0.1768556
# -0.003911017
INV	low	66	2626	47		2241

TRA	high	1	2287	1		2691
# 0.5%     2.5%      50%    97.5%    99.5% 
# -2.2e-03 -1.5e-03  9.9e-05  1.9e-03  2.7e-03 
# 0.5609462
# 0.0001315294
TRA	low	1	2691	1		2287

INV_big high    10      2278    13              2679
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00563 -0.00430 -0.00041  0.00361  0.00502 
# 0.4167848
# -0.0003927225
INV_big low     13      2679    10              2278

y1 = 66
n1 = y1+2626
y2 = 47
n2 = y2+2241
psi_1 = rbeta(10000000, y1+1, (n1-y1)+1)  
psi_2 = rbeta(10000000, y2+1, (n2-y2)+1)
diff_psi1_psi2 = psi_1 - psi_2
diff_quantiles = quantile(diff_psi1_psi2,c(0.005,0.025,0.5,0.975,0.995)) #,na.rm=TRUE)
print(diff_quantiles,digits=2)
print(mean(psi_1 > psi_2))
print(mean(diff_psi1_psi2))

bitmap("fig_1_SVs_on_DEs_new.tiff", height = 8, width = 12, units = 'in', res=1600)
ggplot() + theme_bw() +
  aes(diff_psi1_psi2) +
  geom_density(aes(y = ..scaled..)) +
  geom_vline(aes(xintercept=diff_quantiles[2]), col="blue") +
  geom_vline(aes(xintercept=diff_quantiles[3]), col="red",linetype=2) +
  geom_vline(aes(xintercept=diff_quantiles[4]), col="blue") +
  labs(x = "psi_1 - psi_2 (Differential gene expression occurs with higher chance than no differential expression)" , y="p(psi_1 - psi_2 | y, n)",
       title = "Posterior Probability Distribution of Differential-expression or no change in expression when SV occurs")
dev.off()


#Significant?
#Na Na ! wrong question!!!!
#You should say, based on the observed data, there's a 80% chance there's a 
#higher chance of seeing SVs when there is differential expressions than no DE, or that there’s 
#a 95% chance (Credible interval) the difference between Gene-expression change in SVs and 
#when there is no expression change then prevalence of no-SV falls in the range (-0.0016,0.0086).

