setwd("/scr/k61san/nowicklab/SV-detection/BLAT/compare_sim/runs/Manuscript")

library(ggplot2)
library(gridExtra)


blasr_Del_S<-read.table("./extract/check_extract_single_del_score_top_blasr_all_synt", header=F)
bwa_Del_S<-read.table("./extract/check_extract_single_del_score_top_bwa_all_synt", header=F)
lastz_Del_S<-read.table("./extract/check_extract_single_del_score_top_lastz_all_synt", header=F)

# par(mfrow=c(3,1)) 
# 
# hist(blasr_Del_S$V13)
# hist(bwa_Del_S$V13)
# hist(lastz_Del_S$V13)
# 
# par(mfrow=c(3,1)) 
# h1 <- hist(blasr_Del_S[,13], breaks=50, plot=FALSE)
# cuts <- cut(h1$breaks, c(-Inf,0.5,1,1.5,2,2.5,3,Inf))
# plot(h1, col=cuts, main = "Histogram of Normalised scores with BLASR")
# 
# h2 <- hist(bwa_Del_S[,13], breaks=50, plot=FALSE)
# plot(h2, col=cuts, main = "Histogram of Normalised scores with BWA")
# 
# h3 <- hist(lastz_Del_S[,13], breaks=50, plot=FALSE)
# plot(h3, col=cuts, main = "Histogram of Normalised scores with LASTZ")
# 
# 
# m1 <- ggplot(blasr_Del_S, aes(x = blasr_Del_S[,13]))
# #m1 + geom_histogram(breaks=seq(-1, 4, by =0.5), col="black", aes(fill=..count..))
# m1 + geom_histogram(breaks=seq(-1, 4, by =0.5), col="black", aes(y=..count../sum(..count..),fill = factor(blasr_Del_S[,14])))
# m1 <- ggplot(blasr_Del_S, aes(x = blasr_Del_S[,13])) + 
#   geom_histogram(breaks=seq(-1, 4, by =0.5), col="black", aes(y=..count../sum(..count..),fill = factor(blasr_Del_S[,14])))
# 
# m2 <- ggplot(bwa_Del_S, aes(x = bwa_Del_S[,13]))
# m2 + geom_histogram(breaks=seq(-1, 4, by =0.5), col="black", aes(y=..count../sum(..count..),fill = factor(bwa_Del_S[,14])))
# 
# m3 <- ggplot(lastz_Del_S, aes(x = lastz_Del_S[,13]))
# m3 + geom_histogram(breaks=seq(-1, 4, by =0.5), col="black", aes(y=..count../sum(..count..),fill = factor(lastz_Del_S[,14])))



par(mfrow=c(3,1)) 
m1 <- ggplot(blasr_Del_S, aes(x = blasr_Del_S[,13])) + 
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(blasr_Del_S[,14])))
m2 <- ggplot(bwa_Del_S, aes(x = bwa_Del_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(bwa_Del_S[,14])))
m3 <- ggplot(lastz_Del_S, aes(x = lastz_Del_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(lastz_Del_S[,14])))


grid.arrange(m1,m2,m3, ncol = 1)


blasr_Ins_S<-read.table("./extract/check_extract_single_ins_score_top_blasr_all_synt", header=F)
bwa_Ins_S<-read.table("./extract/check_extract_single_ins_score_top_bwa_all_synt", header=F)
lastz_Ins_S<-read.table("./extract/check_extract_single_ins_score_top_lastz_all_synt", header=F)

I1 <- ggplot(blasr_Ins_S, aes(x = blasr_Ins_S[,13])) + 
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(blasr_Ins_S[,14])))
I2 <- ggplot(bwa_Ins_S, aes(x = bwa_Ins_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(bwa_Ins_S[,14])))
I3 <- ggplot(lastz_Ins_S, aes(x = lastz_Ins_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(lastz_Ins_S[,14])))

grid.arrange(I1,I2,I3, ncol = 1)


blasr_Inv_S<-read.table("./extract/check_extract_single_inv_score_top_blasr_all_synt", header=F)
bwa_Inv_S<-read.table("./extract/check_extract_single_inv_score_top_bwa_all_synt", header=F)
lastz_Inv_S<-read.table("./extract/check_extract_single_inv_score_top_lastz_all_synt", header=F)

Iv1 <- ggplot(blasr_Inv_S, aes(x = blasr_Inv_S[,13])) + 
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(blasr_Inv_S[,14])))
Iv2 <- ggplot(bwa_Inv_S, aes(x = bwa_Inv_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(bwa_Inv_S[,14])))
Iv3 <- ggplot(lastz_Inv_S, aes(x = lastz_Inv_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(lastz_Inv_S[,14])))

grid.arrange(Iv1,Iv2,Iv3, ncol = 1)


blasr_Dup_S<-read.table("./extract/check_extract_single_dup_score_top_blasr_all_synt", header=F)
bwa_Dup_S<-read.table("./extract/check_extract_single_dup_score_top_bwa_all_synt", header=F)
lastz_Dup_S<-read.table("./extract/check_extract_single_dup_score_top_lastz_all_synt", header=F)

Dp1 <- ggplot(blasr_Dup_S, aes(x = blasr_Dup_S[,13])) + 
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(blasr_Dup_S[,14])))
Dp2 <- ggplot(bwa_Dup_S, aes(x = bwa_Dup_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(bwa_Dup_S[,14])))
Dp3 <- ggplot(lastz_Dup_S, aes(x = lastz_Dup_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(lastz_Dup_S[,14])))

grid.arrange(Dp1,Dp2,Dp3, ncol = 1)


blasr_Trans_S<-read.table("./extract/check_extract_single_trans_score_top_blasr_all_synt", header=F)
bwa_Trans_S<-read.table("./extract/check_extract_single_trans_score_top_bwa_all_synt", header=F)
lastz_Trans_S<-read.table("./extract/check_extract_single_trans_score_top_lastz_all_synt", header=F)

T1 <- ggplot(blasr_Trans_S, aes(x = blasr_Trans_S[,13])) + 
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(blasr_Trans_S[,14])))
T2 <- ggplot(bwa_Trans_S, aes(x = bwa_Trans_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(bwa_Trans_S[,14])))
T3 <- ggplot(lastz_Trans_S, aes(x = lastz_Trans_S[,13])) +
  geom_histogram(breaks=seq(-0.5, 3.5, by =0.5), col="black", aes(y=..count..,fill = factor(lastz_Trans_S[,14])))

grid.arrange(T1,T2,T3, ncol = 1)


