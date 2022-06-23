setwd("C:/Users/rohit/Desktop/GenomePaper/SVs/")

library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(reshape2)
library(colorspace)
# library(cowplot)


dfm1<-read.table("Usgn_distances_SVs_to_CDS.txt", header=F, sep = "\t")

mf2 <- split(dfm1, dfm1$V2)

# ggplot(mf2$DEL, aes(mf2$DEL[,3])) + 
#   geom_density(aes(fill=factor(mf2$DEL[,1])), alpha=0.15) + 
#   ##scale_fill_manual( values = c("Positive"="red","Nonpositive"="blue")) +
#   #stat_density(adjust = 2, alpha=0.4) +
#   scale_y_continuous(limits=c("0","5.6e-4"), labels=c(0,1.5e-4,3e-4,4.5e-4),
#                      breaks=c(0,1.5e-4,3e-4,4.5e-4)) + 
#   scale_x_continuous(limits=c("-10000","10000"), labels=c(-10000,-5000,0,5000,10000),
#                      breaks=c(-10000,-5000,0,5000,10000)) + 
#   ggtitle("Deletions") +
#   ylab("Distribution of distances") + xlab("Distance of SV to feature") +
#   theme_bw() + theme(legend.title=element_blank())+
#   annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
#            fill="blue", colour="brown", alpha=0.05)



g1 <- ggplot(mf2$DEL, aes(mf2$DEL[,3])) + 
  geom_density(aes(fill=factor(mf2$DEL[,1])), alpha=0.15) + 
  ##scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1.3e-3) + xlim(-2000,7000) + ggtitle("Deletions (Overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g2 <- ggplot(mf2$DUP, aes(mf2$DUP[,3])) + 
  geom_density(aes(fill=factor(mf2$DUP[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1.1e-3) + xlim(-8000,30000) + ggtitle("Duplications (Overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


g3 <- ggplot(mf2$INS, aes(mf2$INS[,3])) + 
  geom_density(aes(fill=factor(mf2$INS[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,0.75e-3) + xlim(-3500,28000) + ggtitle("Insertions (Overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g4 <- ggplot(mf2$INV, aes(mf2$INV[,3])) + 
  geom_density(aes(fill=factor(mf2$INV[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,0.51e-3) + xlim(-8000,25000) + ggtitle("Inversions (Overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g5 <- ggplot(mf2$TRA, aes(mf2$TRA[,3])) + 
  geom_density(aes(fill=factor(mf2$TRA[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1.5e-5) + xlim(-50000,100000) + ggtitle("Translocations (Overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


grid.arrange(g1,g2,g3,g4,ncol=2)


dat1_1 <- subset(dfm1,V1=="Positive")

g_P_1 <- ggplot(dat1_1, aes(dat1_1[,3])) + 
  geom_density(aes(fill=factor(dat1_1[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,3e-4) + xlim(-5000,25000) + ggtitle("Influence of SVs on Positively selected regions (overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())


dat1_2 <- subset(dfm1,V1=="Nonpositive")

g_N_1 <- ggplot(dat1_2, aes(dat1_2[,3])) + 
  geom_density(aes(fill=factor(dat1_2[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1e-4) + xlim(-2000,35000) + ggtitle("Influence of SVs on Nonpositive regions (overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())


grid.arrange(g_P_1,g_N_1,ncol=1)

dfm2<-read.table("Usgn_noovl_distances_SVs_to_CDS.txt", header=F, sep = "\t")

mf3 <- split(dfm2, dfm2$V2)


g6 <- ggplot(mf3$DEL, aes(mf3$DEL[,3])) + 
  geom_density(aes(fill=factor(mf3$DEL[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,8.75e-5) + xlim(-10000,40000) + ggtitle("Deletions (Without overlaps)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g7 <- ggplot(mf3$DUP, aes(mf3$DUP[,3])) + 
  geom_density(aes(fill=factor(mf3$DUP[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,7.5e-5) + xlim(-10000,40000) + ggtitle("Duplications (Without overlaps)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


g8 <- ggplot(mf3$INS, aes(mf3$INS[,3])) + 
  geom_density(aes(fill=factor(mf3$INS[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1e-4) + xlim(-10000,50000) + ggtitle("Insertions (Without overlaps)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

g9 <- ggplot(mf3$INV, aes(mf3$INV[,3])) + 
  geom_density(aes(fill=factor(mf3$INV[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,9e-5) + xlim(-10000,60000) + ggtitle("Inversions (Without overlaps)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


g10 <- ggplot(mf3$TRA, aes(mf3$TRA[,3])) + 
  geom_density(aes(fill=factor(mf3$TRA[,1])), alpha=0.15) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,2.6e-5) + xlim(-10000,75000) + ggtitle("Translocations (Without overlaps)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())



grid.arrange(g6,g7,g8,g9,ncol=2)


dat2_1 <- subset(dfm2,V1=="Positive")

g_P_2 <- ggplot(dat2_1, aes(dat2_1[,3])) + 
  geom_density(aes(fill=factor(dat2_1[,1])), alpha=0.15) + 
  scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,7.5e-5) + xlim(-10000,100000) + ggtitle("Influence of non-overlapping SVs on Positively selected regions") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


dat2_2 <- subset(dfm2,V1=="Nonpositive")

g_N_2 <- ggplot(dat2_2, aes(dat2_2[,3])) + 
  geom_density(aes(fill=factor(dat2_2[,1])), alpha=0.15) + 
  stat_density(adjust = 2, alpha=0.4) + 
  #scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  ylim(0,9e-5) + xlim(-10000,50000) + ggtitle("Influence of non-overlapping SVs on Nonpositive regions") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)


grid.arrange(g_P_2,g_N_2,ncol=1)


grid.arrange(g1,g2,g3,g4,g5,ncol=2)
grid.arrange(g_P_1,g_N_1,ncol=1)
grid.arrange(g6,g7,g8,g9,g10,ncol=2)
grid.arrange(g_P_2,g_N_2,ncol=1)

grid.arrange(g_P_1,g_N_1,g_P_2,g_N_2,ncol=2)

svg("proximity_withovl_unsigned.svg", height = 8, width = 12, pointsize=12)
grid.arrange(g1,g2,g3,g4,ncol=2,  
             sub=textGrob("A) Absolute distances of overlapping and proximity of SVs to coding regions", vjust=0, 
                          gp = gpar(fontsize=14,   fontfamily="Times New Roman")))
dev.off()


grid.arrange(g_P_1,g_N_1,ncol=1)

svg("proximity_noovl_unsigned.svg", height = 8, width = 12, pointsize=12)
grid.arrange(g6,g7,g8,g9,ncol=2,  
             sub=textGrob("B) Absolute distances of proximity of SVs to coding regions without considering overlaps", vjust=0, 
                          gp = gpar(fontsize=14,   fontfamily="Times New Roman")))
dev.off()



grid.arrange(g_P_2,g_N_2,ncol=1)


grid.arrange(g5,g10,ncol=1, 
             sub=textGrob("Absolute distances of proximity of Translocations to coding regions with and without overlaps",  
                          gp =gpar(fontsize=14, fontfamily="Times New Roman"), vjust=0))


grid.arrange(g_P_1,g_N_1,g_P_2,g_N_2,ncol=2,  
             sub=textGrob("Absolute distances of proximity of coding regions to the closest SV categorised by the type of selection", vjust=0, 
                          gp = gpar(fontsize=14,   fontfamily="Times New Roman")))




G1 <- ggplot(mf2$DEL, aes(mf2$DEL[,3])) + 
  geom_density(aes(fill=factor(mf2$DEL[,1])), alpha=0.15) + 
  ##scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ylim(0,1.3e-3) + xlim(-2000,7000) + ggtitle("Deletions (Overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())+
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, 
           fill="green", colour="brown", alpha=0.05)

G2 <- ggplot(mf2$DEL, aes(x=mf2$DEL[,3], fill=mf2$DEL[,1])) + 
  geom_histogram(alpha=0.15) + 
  ##scale_fill_manual( values = c("Positive"="#00BA38","Nonpositive"="#F8766D")) +
  #stat_density(adjust = 2, alpha=0.4) +
  ggtitle("Deletions (Overlaps accounted)") +
  ylab("Distribution of distances") + xlab("Unsigned distance of SV to feature") +
  theme_bw() + theme(legend.title=element_blank())


grid.arrange(G1,G2,ncol=2)
