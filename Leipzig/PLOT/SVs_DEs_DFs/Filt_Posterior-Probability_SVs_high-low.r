
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

complex low     3       2689    1               2287
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00344 -0.00261 -0.00059  0.00131  0.00213 
# > print(mean(psi_1 > psi_2))
# [1] 0.2421998
# > print(mean(diff_psi1_psi2))
# [1] -0.000611341
DEL     low     1110    1582    881             1407
# 0.5%     2.5%      50%    97.5%    99.5% 
# -8.6e-03 -2.4e-05  2.7e-02  5.4e-02  6.3e-02 
# > print(mean(psi_1 > psi_2))
# [1] 0.9748965
# > print(mean(diff_psi1_psi2))
# [1] 0.02725016
DUP     low     110     2582    68              2220
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0025  0.0008  0.0111  0.0213  0.0246 
# > print(mean(psi_1 > psi_2))
# [1] 0.9825526
# > print(mean(diff_psi1_psi2))
# [1] 0.0110737
IDP     low     1       2691    1               2287
# 0.5%     2.5%      50%    97.5%    99.5% 
# -2.7e-03 -1.9e-03 -9.9e-05  1.5e-03  2.2e-03 
# > print(mean(psi_1 > psi_2))
# [1] 0.4392183
# > print(mean(diff_psi1_psi2))
# [1] -0.0001313482
INS     low     619     2073    507             1781
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0223 -0.0150  0.0083  0.0316  0.0389 
# > print(mean(psi_1 > psi_2))
# [1] 0.7582395
# > print(mean(diff_psi1_psi2))
# [1] 0.008311675
INVDUP  low     0       2692    0               2288
# INV     low     66      2626    47              2241
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0071 -0.0044  0.0039  0.0122  0.0148 
# > print(mean(psi_1 > psi_2))
# [1] 0.8230138
# > print(mean(diff_psi1_psi2))
# [1] 0.003909063
TRA     low     1       2691    1               2287


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

###

complex	high	1	2287	0		103
complex	low	0	129	1		2261
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00185 -0.00097  0.00446  0.02714  0.03908 
# > print(mean(psi_1 > psi_2))
# [1] 0.8942138
# > print(mean(diff_psi1_psi2))
# [1] 0.006751979

DEL	high	881	1407	44		59
DEL	low	53	76	872		1390
# 0.5%   2.5%    50%  97.5%  99.5% 
# -0.084 -0.058  0.026  0.114  0.142 
# > print(mean(psi_1 > psi_2))
# [1] 0.724255
# > print(mean(diff_psi1_psi2))
# [1] 0.02660821

DUP	high	68	2220	8		95
DUP	low	9	120	67		2195
# 0.5%    2.5%     50%   97.5%   99.5% 
# -0.0021  0.0067  0.0442  0.0977  0.1179 
# > print(mean(psi_1 > psi_2))
# [1] 0.9923201
# > print(mean(diff_psi1_psi2))
# [1] 0.04629153

IDP	high	1	2287	0		103
IDP	low	0	129	1		2261
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00186 -0.00097  0.00446  0.02710  0.03909 
# > print(mean(psi_1 > psi_2))
# [1] 0.8941684
# > print(mean(diff_psi1_psi2))
# [1] 0.006748056

INS	high	507	1781	24		79
INS	low	33	96	498		1764
# 0.5%   2.5%    50%  97.5%  99.5% 
# -0.055 -0.034  0.038  0.119  0.146 
# > print(mean(psi_1 > psi_2))
# [1] 0.8408964
# > print(mean(diff_psi1_psi2))
# [1] 0.039141

INVDUP	high	0	2288	0		103
INVDUP	low	0	129	0		2262

INV	high	47	2241	8		95
INV	low	8	121	47		2215
# 0.5%   2.5%    50%  97.5%  99.5% 
# 0.0025 0.0104 0.0454 0.0968 0.1164 
# > print(mean(psi_1 > psi_2))
# [1] 0.9973422
# > print(mean(diff_psi1_psi2))
# [1] 0.04749819

TRA	high	1	2287	0		103
TRA	low	0	129	1		2261
# 0.5%     2.5%      50%    97.5%    99.5% 
# -0.00186 -0.00097  0.00446  0.02712  0.03905 
# > print(mean(psi_1 > psi_2))
# [1] 0.8940829
# > print(mean(diff_psi1_psi2))
# [1] 0.006751024

y1 = 0
n1 = y1+129
y2 = 1
n2 = y2+2261
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


