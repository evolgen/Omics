
setwd("/homes/biertank/rohit/Downloads/scripts/SVs_DEs_DFs")

# Libraries you need
library(ggplot2)

# Loading Data
grfTable <- read.table("/homes/biertank/rohit/Downloads/scripts/SVs_DEs_DFs/in_Table",sep = "\t",header = TRUE)

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


#### SVs and no-SVs with DE and no-DE
y1 = 916
y2 = 7001
n1 = 3775+916
n2 = 24210+7001

# 0.5%   2.5%    50%  97.5%  99.5% 
# -0.045 -0.041 -0.029 -0.017 -0.013 
# 2.7e-06
# -0.02893198



#### SVs and no-SVs with high-DE and no-high-DE
y1 = 1187
n1 = 1187+8326
y2 = 3467
n2 = 3467+28575

# 0.5%  2.5%   50% 97.5% 99.5% 
# 0.0069 0.0092 0.0166 0.0241 0.0266 
# 
# 0.9999953
# 
# 0.01663018


#### SVs and no-SVs with low-DE and no-low-DE
y1 = 1609
n1 = 9519
y2 = 3903
n2 = 32097


# 0.5%  2.5%   50% 97.5% 99.5% 
# 0.037 0.039 0.047 0.056 0.059 
# 
# 1
# 
# 0.04747626


#### SVs and no-SVs with high-DE and low-DE
y1 = 1187
n1 = 1187+2604
y2 = 3467
n2 = 3467+3903

# 0.5%  2.5%   50% 97.5% 99.5% 
# -0.18 -0.18 -0.16 -0.14 -0.13 
# 
# 0
# 
# -0.1572178


# Deletions and SV
y1 = 2259
y2 = 7666
n1 = 2259+5078
n2 = 7666+24764
# 0.5%  2.5%   50% 97.5% 99.5% 
# 0.056 0.060 0.072 0.083 0.087 
# 1
# 0.071543


# Insertions and SV
y1 = 506
y2 = 8103
n1 = 506+1034
n2 = 8103+25660
# 0.5%  2.5%   50% 97.5% 99.5% 
# 0.058 0.065 0.089 0.113 0.121 
# 1
# 0.08878621


# Inversions and SV
y1 = 459
y2 = 8125
n1 = 459+1292
n2 = 8125+25698
# 0.5%  2.5%   50% 97.5% 99.5% 
# -0.0050  0.0013  0.0221  0.0435  0.0504 
# 0.9816453
# 0.022172


# Duplications and SV
y1 = 587
y2 = 7520
n1 = 587+1477
n2 = 7520+24221
# 0.5%  2.5%   50% 97.5% 99.5% 
# 0.022 0.028 0.048 0.068 0.074 
# 0.9999994
# 0.04767456


psi_1 = rbeta(10000000, y1+1, (n1-y1)+1)  
psi_2 = rbeta(10000000, y2+1, (n2-y2)+1)
diff_psi1_psi2 = psi_1 - psi_2

diff_quantiles = quantile(diff_psi1_psi2,c(0.005,0.025,0.5,0.975,0.995)) #,na.rm=TRUE)
print(diff_quantiles,digits=2)

print("For the occurance of SV, the probability of higher chance in DE than no DE occurence:")
print(mean(psi_1 > psi_2))

print("Mean difference psi_1-psi_2")
print(mean(diff_psi1_psi2))

bitmap("fig_1_SVs_on_DEs.tiff", height = 8, width = 12, units = 'in', res=1600)
ggplot() + theme_bw() +
  aes(diff_psi1_psi2) +
  geom_density(aes(y = ..scaled..)) +
  geom_vline(aes(xintercept=diff_quantiles[2]), col="blue") +
  geom_vline(aes(xintercept=diff_quantiles[3]), col="red",linetype=2) +
  geom_vline(aes(xintercept=diff_quantiles[4]), col="blue") +
  labs(x = "psi_1 - psi_2 (Differential gene expression occurs with higher chance than no differential expression)" , y="p(psi_1 - psi_2 | y, n)",
       title = "Posterior Probability Distribution of Differential-expression or no change in expression when SV occurs")
dev.off()

y1 = 354
y2 = 3940
n1 = 354+1723
n2 = 3940+27801

# SIMULATION based on Beta distribution
# 1000000 is the number of random points and should be enough! Higher number gives higher accuracy (more simulations)
psi_1 = rbeta(10000000, y1+1, (n1-y1)+1)  
psi_2 = rbeta(10000000, y2+1, (n2-y2)+1)
diff_psi1_psi2 = psi_1 - psi_2  # simulated difference

# SIMULATION based on Dirichlet distribution
# This is one level zoom out from the Domain change aand saying given a 
# Change in one Gene does it more human specific or not or any other question
# Sanity check: If we only look at one variable or others, it means we are using Beta
#Dtheta1=rdirichlet(I,c(y1+1, (n1-y1)+1))
#Dtheta2=rdirichlet(I,c(y2+1, (n2-y2)+1))
#Ddiff = Dtheta1-Dtheta2  # simulated diffs

# OUTPUT
#This numbers are showing the posterior quantiles for the difference. 
#The posterior median (50% quantile) difference is 0.116. The other numbers are 
#quantiles which determine 95% and 99% posterior intervals, which span (0.025,0.975) 
#and (0.005,0.995) respectively. The posterior central 95% quantile is thus (-0.060,0.355). 

diff_quantiles = quantile(diff_psi1_psi2,c(0.005,0.025,0.5,0.975,0.995)) #,na.rm=TRUE)
print(diff_quantiles,digits=2)

#The number will print is the posterior probability that Gene-expression change occures at a 
#a higher percentage of due to the presence of SV compared to no-SV, which is 1. This is also the value 
#you would get from integrating the posterior density from 0 to \infty.

print("For the occurance of SV, the probability of higher chance in DE than no DE occurence:")
print(mean(psi_1 > psi_2))

#The final number is the mean difference between psi_1 and psi_2, which is 
#the unbiased Bayesian estimator of the difference.

print("Mean difference psi_1-psi_2")
print(mean(diff_psi1_psi2))

# VISUALIZATION
# generates the graph of the Posterior Probability Distribution. The blue vertical 
# lines are drawn at the boundaries of the central 95% region of the posterior 
# difference p(\psi_1 - \psi_2 | n, y).

ggplot() + theme_bw() +
  aes(diff_psi1_psi2) +
  geom_density(aes(y = ..scaled..)) +
  geom_vline(aes(xintercept=diff_quantiles[2]), col="blue") +
  geom_vline(aes(xintercept=diff_quantiles[3]), col="red",linetype=2) +
  geom_vline(aes(xintercept=diff_quantiles[4]), col="blue") +
  labs(x = "psi_1 - psi_2 (Differential gene expression occurs more than no differential expression)" , y="p(psi_1 - psi_2 | y, n)",
       title = "Posterior Probability Distribution of Differential-expression or no DE when SV occurs")

#Significant?
#Na Na ! wrong question!!!!
#You should say, based on the observed data, there's a 80% chance there's a 
#higher chance of seeing SVs when there is differential expressions than no DE, or that there’s 
#a 95% chance (Credible interval) the difference between Gene-expression change in SVs and 
#when there is no expression change then prevalence of no-SV falls in the range (-0.0016,0.0086).


