
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

y1=0
x1=4
n1=x1+y1
y2=3
x2=34
n2=x2+y2

psi_1 = rbeta(10000000, y1+1, (n1-y1)+1)  
psi_2 = rbeta(10000000, y2+1, (n2-y2)+1)
diff_psi1_psi2 = psi_1 - psi_2
diff_quantiles = quantile(diff_psi1_psi2,c(0.005,0.025,0.5,0.975,0.995)) #,na.rm=TRUE)
print(diff_quantiles,digits=2)
print(mean(psi_1 > psi_2))
print(mean(diff_psi1_psi2))

complex positive        0       772     0               6014
complex neutral 0       6253    0               546
DEL     positive        220     552     2625            3389
DEL     neutral 2723    3530    156             390
DUP     positive        27      745     238             5776
DUP     neutral 250     6003    19              527
IDP     positive        0       772     2               6012
IDP     neutral 2       6251    0               546
INS     positive        130     642     1660            4354
INS     neutral 1711    4542    94              452
INVDUP  positive        0       772     0               6014
INVDUP  neutral 0       6253    0               546
INV     positive        18      754     162             5852
INV     neutral 170     6083    12              534
TRA     positive        1       771     2               6012
TRA     neutral 2       6251    1               545




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


