setwd("C:/Users/rohit/Desktop/GenomePaper/SVs/")

library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(reshape2)
library(colorspace)
#library(cowplot)


dfm1<-read.table("distances_SVs_to_CDS.txt", header=F, sep = "\t")

# dat2 <- dfm1[,c(1,2,3)]
# 
# dfm2 <- melt(dat2, id.vars=c(1,2),measure.vars = c(3))


#bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/Repeatmasker_species_tetrapods_2.tiff", 
#       height = 8, width = 12, units = 'in', res=1200)

pal <- diverge_hcl(7)

# 
# ggplot(dfm1, aes(dfm1$V1,dfm1$V2)) +
#   geom_point(stat = "identity")
# 
# pairs(~dfm1$V2+dfm1$V2,data=dfm1, 
#       main="Simple Scatterplot Matrix")
# 
# plot(dfm1$V2, dfm1$V3, pch=21, 
#     bg=c("red","green3")[unclass(dfm1$v1)], 
#     main="check)

# ggplot(dfm1, aes(x = dfm1$V2) ) + geom_density(color=factor(dfm1$V1)) 
# 
# mf <- split(dfm1, dfm1$V1)
# 
# m1 <- ggplot(mf$Positive, aes(y = log10(mf$Positive[,3]), x = seq(1, length(mf$Positive[,3])))) + 
#   geom_point(aes(color=factor(mf$Positive[,2]), shape=factor(mf$Positive[,2])))
# 
# 
# m2 <- ggplot(mf$Neutral, aes(y = log10(mf$Neutral[,3]), x = seq(1, length(mf$Neutral[,3])))) + 
#   geom_point(aes(color=factor(mf$Neutral[,2])))
# 
# grid.arrange(m1,m2, ncol=1)
# 
# m1 <- ggplot(mf$Positive, aes(log10(mf$Positive[,3]))) + 
#   geom_density(aes(fill=factor(mf$Positive[,2]), alpha=0.01))
# 
# m2 <- ggplot(mf$Neutral, aes(mf$Neutral[,3])) + 
#   geom_density(aes(fill=factor(mf$Neutral[,2]), alpha=0.01)) + xlim(-5,5)
# 
# grid.arrange(m1,m2, ncol=1)
# 
# ggplot(mf$Positive, aes(mf$Positive[,3])) + 
#   geom_density(aes(fill=factor(mf$Positive[,2]), alpha=2)) + 
#   stat_density(adjust = 2, alpha=0.4) +
#   ylim(0,6.7e-5) + xlim(-40000,40000)
#     scale_y_continuous(trans="log10")
# 
# 
# ggplot(mf$Neutral, aes(mf$Neutral[,3])) + 
#   geom_density(aes(fill=factor(mf$Neutral[,2]), alpha=0.01)) + 
#   stat_density(adjust = 2, alpha=0.5) +
#   stat_density(adjust = 2, alpha=0.1) +
#   ylim(0,2.5e-4) + xlim(-40000,50000)
# 
# scale_y_continuous(trans="log10", name="density") 

mf2 <- split(dfm1, dfm1$V2)


g1 <- ggplot(mf2$DEL, aes(mf2$DEL[,3])) + 
  geom_density(aes(fill=factor(mf2$DEL[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,5.6e-4) + xlim(-10000,10000) + ggtitle("Deletions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g2 <- ggplot(mf2$DUP, aes(mf2$DUP[,3])) + 
  geom_density(aes(fill=factor(mf2$DUP[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,3.2e-4) + xlim(-7000,20000) + ggtitle("Duplications") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


g3 <- ggplot(mf2$INS, aes(mf2$INS[,3])) + 
  geom_density(aes(fill=factor(mf2$INS[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,2e-4) + ggtitle("Insertions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+ 
  scale_x_continuous(limits=c(-30000,30000),
  breaks=c(-25000,-20000,-15000,-10000,-5000,-1000,1000,5000,10000,15000,20000,25000)) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05) +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank())

g4 <- ggplot(mf2$INV, aes(mf2$INV[,3])) + 
  geom_density(aes(fill=factor(mf2$INV[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1.9e-4) + ggtitle("Inversions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  scale_x_continuous(limits=c(-20000,20000),
     breaks=c(-10000,-5000,-1000,1000,5000,10000,15000,20000)) +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


g5 <- ggplot(mf2$TRA, aes(mf2$TRA[,3])) + 
  geom_density(aes(fill=factor(mf2$TRA[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1.75e-5) + xlim(-100000,100000) + ggtitle("Translocations (with overlaps)") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)



grid.arrange(g1,g2,g3,g4,g5,ncol=2)


dat1_1 <- subset(dfm1,V1=="Positive")

g_P_1 <- ggplot(dat1_1, aes(dat1_1[,3])) + 
  geom_density(aes(fill=factor(dat1_1[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,2.2e-4) + xlim(-18000,18000) + ggtitle("Influence of SVs on Positively selected regions (overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
   theme_bw() + theme(legend.title=element_blank()) #+
  # annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
  #          fill="green", colour="brown", alpha=0.05)


dat1_2 <- subset(dfm1,V1=="Neutral")

g_N_1 <- ggplot(dat1_2, aes(dat1_2[,3])) + 
  geom_density(aes(fill=factor(dat1_2[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1.4e-3) + xlim(-2000,5000) + ggtitle("Influence of SVs on Neutral regions (overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())


grid.arrange(g_P_1,g_N_1,ncol=1)

dfm2<-read.table("novl_distances_SVs_to_CDS.txt", header=F, sep = "\t")

mf3 <- split(dfm2, dfm2$V2)


g6 <- ggplot(mf3$DEL, aes(mf3$DEL[,3])) + 
  geom_density(aes(fill=factor(mf3$DEL[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,9e-5) + xlim(-35000,55000) + ggtitle("Deletions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g7 <- ggplot(mf3$DUP, aes(mf3$DUP[,3])) + 
  geom_density(aes(fill=factor(mf3$DUP[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,4e-5) + xlim(-50000,55000) + ggtitle("Duplications") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


g8 <- ggplot(mf3$INS, aes(mf3$INS[,3])) + 
  geom_density(aes(fill=factor(mf3$INS[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,6e-5) + xlim(-50000,50000) + ggtitle("Insertions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g9 <- ggplot(mf3$INV, aes(mf3$INV[,3])) + 
  geom_density(aes(fill=factor(mf3$INV[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,5.5e-5) + xlim(-65000,65000) + ggtitle("Inversions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


g10 <- ggplot(mf3$TRA, aes(mf3$TRA[,3])) + 
  geom_density(aes(fill=factor(mf3$TRA[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1e-5) + xlim(-200000,200000) + ggtitle("Translocations (No overlaps)") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


grid.arrange(g6,g7,g8,g9,g10,ncol=2)


dat2_1 <- subset(dfm2,V1=="Positive")

g_P_2 <- ggplot(dat2_1, aes(dat2_1[,3])) + 
  geom_density(aes(fill=factor(dat2_1[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  stat_density(adjust = 2, alpha=0.4) +
  ylim(0,5.5e-5) + xlim(-30000,75000) + ggtitle("Influence of non-overlapping SVs on Positively selected regions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.15) #+


dat2_2 <- subset(dfm2,V1=="Neutral")

g_N_2 <- ggplot(dat2_2, aes(dat2_2[,3])) + 
  geom_density(aes(fill=factor(dat2_2[,1])), alpha=0.15) + 
  stat_density(adjust = 2, alpha=0.4) + 
  scale_fill_manual( values = c("Positive"="red","Neutral"="blue")) +
  ylim(0,5.7e-5) + xlim(-25000,50000) + ggtitle("Influence of non-overlapping SVs on Neutral regions") +
  ylab("Distribution of distances") + xlab("Distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


grid.arrange(g_P_2,g_N_2,ncol=1)


grid.arrange(g1,g2,g3,g4,ncol=2,  
  sub=textGrob("A) Signed distances of overlapping and proximity of SVs to coding regions", vjust=0, 
  gp = gpar(fontsize=14,   fontfamily="Times New Roman")))


grid.arrange(g_P_1,g_N_1,ncol=1)

grid.arrange(g6,g7,g8,g9,ncol=2,  
             sub=textGrob("B) Signed distances of proximity of SVs to coding regions without considering overlaps", vjust=0, 
                          gp = gpar(fontsize=14,   fontfamily="Times New Roman")))

grid.arrange(g_P_2,g_N_2,ncol=1)


grid.arrange(g5,g10,ncol=1, 
  sub=textGrob("Signed distances of proximity of Translocations to coding regions with and without overlaps",  
  gp =gpar(fontsize=14, fontfamily="Times New Roman"), vjust=0))


grid.arrange(g_P_1,g_N_1,g_P_2,g_N_2,ncol=2,  
             sub=textGrob("Signed distances of proximity of coding regions to the closest SV categorised by the type of selection", vjust=0, 
                          gp = gpar(fontsize=14,   fontfamily="Times New Roman")))
             

dev.off()

# grid.arrange(arrangeGrob(g_P_1,g_N_1,g_P_2,g_N_2, ncol=2,
#             sub=textGrob("A (hr)", vjust=0, rot=90, hjust = 3,
#                          gp = gpar(fontsize=20,   fontfamily="Times New Roman")),
#             left=textGrob("B (MPH)", rot=90, hjust = 1,
#                           gp =gpar(fontsize=18, fontfamily="Times New Roman"), vjust=0)))
# 
# 
