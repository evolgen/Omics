setwd("/scr/k61san/nowicklab/SV-detection/BLAT/compare_sim/runs/Manuscript/extract")

library(ggplot2)
library(gridExtra)


blasr_Del_S<-read.table("check_extract_single_del_score_all_blasr_sim.psl", header=F)
bwa_Del_S<-read.table("check_extract_single_del_score_all_bwa_sim.psl", header=F)
lastz_Del_S<-read.table("check_extract_single_del_score_all_lastz.psl", header=F)

m1 <- ggplot(blasr_Del_S, aes(x = log(blasr_Del_S[,23]))) + 
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(blasr_Del_S[,24])))
m2 <- ggplot(bwa_Del_S, aes(x = log(bwa_Del_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(bwa_Del_S[,24])))
m3 <- ggplot(lastz_Del_S, aes(x = log(lastz_Del_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(lastz_Del_S[,24])))


grid.arrange(m1,m2,m3, ncol = 1)
dev.copy(png, "logscores_distribution_all_single_deletion.png")
dev.off()


blasr_Ins_S<-read.table("check_extract_single_ins_score_all_blasr_sim.psl", header=F)
bwa_Ins_S<-read.table("check_extract_single_ins_score_all_bwa_sim.psl", header=F)
lastz_Ins_S<-read.table("check_extract_single_ins_score_all_lastz.psl", header=F)

I1 <- ggplot(blasr_Ins_S, aes(x = log(blasr_Ins_S[,23]))) + 
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(blasr_Ins_S[,24])))
I2 <- ggplot(bwa_Ins_S, aes(x = log(bwa_Ins_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(bwa_Ins_S[,24])))
I3 <- ggplot(lastz_Ins_S, aes(x = log(lastz_Ins_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(lastz_Ins_S[,24])))

grid.arrange(I1,I2,I3, ncol = 1)
dev.copy(png, "logscores_distribution_all_single_insertion.png")
dev.off()


blasr_Inv_S<-read.table("check_extract_single_inv_score_all_blasr_sim.psl", header=F)
bwa_Inv_S<-read.table("check_extract_single_inv_score_all_bwa_sim.psl", header=F)
lastz_Inv_S<-read.table("check_extract_single_inv_score_all_lastz.psl", header=F)

Iv1 <- ggplot(blasr_Inv_S, aes(x = log(blasr_Inv_S[,23]))) + 
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(blasr_Inv_S[,24])))
Iv2 <- ggplot(bwa_Inv_S, aes(x = log(bwa_Inv_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(bwa_Inv_S[,24])))
Iv3 <- ggplot(lastz_Inv_S, aes(x = log(lastz_Inv_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(lastz_Inv_S[,24])))

grid.arrange(Iv1,Iv2,Iv3, ncol = 1)
dev.copy(png, "logscores_distribution_all_single_inversion.png")
dev.off()


blasr_Dup_S<-read.table("check_extract_single_dup_score_all_blasr_sim.psl", header=F)
bwa_Dup_S<-read.table("check_extract_single_dup_score_all_bwa_sim.psl", header=F)
lastz_Dup_S<-read.table("check_extract_single_dup_score_all_lastz.psl", header=F)

Dp1 <- ggplot(blasr_Dup_S, aes(x = log(blasr_Dup_S[,23]))) + 
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(blasr_Dup_S[,24])))
Dp2 <- ggplot(bwa_Dup_S, aes(x = log(bwa_Dup_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(bwa_Dup_S[,24])))
Dp3 <- ggplot(lastz_Dup_S, aes(x = log(lastz_Dup_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(lastz_Dup_S[,24])))

grid.arrange(Dp1,Dp2,Dp3, ncol = 1)
dev.copy(png, "logscores_distribution_all_single_duplication.png")
dev.off()


blasr_Trans_S<-read.table("check_extract_single_trans_score_all_blasr_sim.psl", header=F)
bwa_Trans_S<-read.table("check_extract_single_trans_score_all_bwa_sim.psl", header=F)
lastz_Trans_S<-read.table("check_extract_single_trans_score_all_lastz.psl", header=F)

T1 <- ggplot(blasr_Trans_S, aes(x = log(blasr_Trans_S[,23]))) + 
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(blasr_Trans_S[,24])))
T2 <- ggplot(bwa_Trans_S, aes(x = log(bwa_Trans_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(bwa_Trans_S[,24])))
T3 <- ggplot(lastz_Trans_S, aes(x = log(lastz_Trans_S[,23]))) +
  geom_histogram(bins=15,col="black", aes(y=..count..,fill = factor(lastz_Trans_S[,24])))

grid.arrange(T1,T2,T3, ncol = 1)
dev.copy(png, "logscores_distribution_all_single_transposition.png")
dev.off()



tiff("/homes/biertank/rohit/Downloads/scripts/simulation_SVs/logscores_distribution_allhits_single_SV.tif", 
     units="in", width=20, height=12.5, res=80)
grid.arrange(m1,m2,m3,I1,I2,I3,Iv1,Iv2,Iv3,T1,T2,T3,Dp1,Dp2,Dp3, ncol = 3)
dev.off()

