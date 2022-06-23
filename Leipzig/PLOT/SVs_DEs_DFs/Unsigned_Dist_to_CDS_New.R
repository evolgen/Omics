setwd("/scr/bloodymary/rohit/Lacerta_viridis/Selection")

library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
#library(reshape2)
library(colorspace)


dfm1<-read.table("unsigned_dist_ALL_Filt", header=F, sep = "\t")

mf2 <- split(dfm1, dfm1$V2)

# ggplot(mf2$DEL, aes(x=mf2$DEL$V3, fill=mf2$DEL$V1)) + geom_histogram(aes(fill=factor(mf2$DEL$V1)), alpha=0.2) 
# ggplot(mf2$DEL, aes(x=mf2$DEL$V3, fill=mf2$DEL$V1)) + geom_histogram(alpha=0.2) 
# ggplot(mf2$DUP, aes(x=mf2$DUP$V3, fill=mf2$DUP$V1)) + geom_area(aes(y = ..count..), stat = "bin", alpha=0.4)
# ggplot(mf2$DEL, aes(mf2$DEL[,3])) + geom_density(aes(fill=factor(mf2$DEL[,1])), alpha=0.15)

svg("/homes/biertank/rohit/Downloads/scripts/SVs_DEs_DFs/Total_unsigned_distances.svg", 
    height = 8, width = 12, pointsize=12)

G0 <- ggplot(dfm1, aes(x=V2, ..count..), fill=dfm1$V1) + scale_y_log10() + 
  geom_bar(aes(fill=dfm1$V1), position = "dodge") +
  ggtitle("Distance of CDS to Deletions") +
  #geom_text(stat='count',aes(label=..count..),position = position_dodge(width = 0.5)) +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  theme_bw() + theme(legend.title=element_blank())
  
dev.off()

G1 <- ggplot(mf2$DEL, aes(x=V3,y=log10(..count..))) + #scale_y_log10() + #scale_y_log10(limits = c(-0.1,NA)) + 
  geom_histogram(data=subset(mf2$DEL, mf2$DEL$V1=='Positive'), 
                 aes(y=..count../sum(..count..), fill="Positive"), alpha=0.3) + 
  geom_histogram(data=subset(mf2$DEL, mf2$DEL$V1=='Non-positive'), 
                 aes(y=..count../sum(..count..), fill="Non-positive"), alpha=0.3) +
  scale_fill_manual(values = c("Positive"="green","Non-positive"="red"), name="Selection", 
                    labels = c("Positive"="Positive","Non-positive"="Non-positive")) +
  # scale_color_manual(values = c("Positive"="green","Non-positive"="red"),name="Selection", alpha=0.3,
  #                    labels = c("Positive"="Positive","Non-positive"="Non-positive"),guide="legend") +
  ggtitle("Distance of CDS to Deletions") +
  ylab("Frequency of distances normalised by counts per bin-length") + xlab("Distance of SV to CDS feature") +
  theme_bw() + #theme(legend.title=element_blank()) + expand_limits(y=1) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)


G2 <- ggplot(mf2$INS, aes(x=V3,y=log10(..count..))) + #scale_y_log10() + #scale_y_log10(limits = c(-0.1,NA)) + 
  geom_histogram(data=subset(mf2$INS, mf2$INS$V1=='Positive'), 
                 aes(y=..count../sum(..count..), fill="Positive"), alpha=0.3) + 
  geom_histogram(data=subset(mf2$INS, mf2$INS$V1=='Non-positive'), 
                 aes(y=..count../sum(..count..), fill="Non-positive"), alpha=0.3) +
  scale_fill_manual(values = c("Positive"="green","Non-positive"="red"), name="Selection", 
                    labels = c("Positive"="Positive","Non-positive"="Non-positive")) +
  # scale_color_manual(values = c("Positive"="green","Non-positive"="red"),name="Selection", alpha=0.3,
  #                    labels = c("Positive"="Positive","Non-positive"="Non-positive"),guide="legend") +
  ggtitle("Distance of CDS to Deletions") +
  ylab("Frequency of distances normalised by counts per bin-length") + xlab("Distance of SV to CDS feature") +
  theme_bw() + #theme(legend.title=element_blank()) + expand_limits(y=1) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)

svg("/homes/biertank/rohit/Downloads/scripts/SVs_DEs_DFs/Indels_unsigned_distances.svg", 
    height = 8, width = 12, pointsize=12)

grid.arrange(G1,G2, ncol=2)
dev.off()




G3 <- ggplot(mf2$DUP, aes(x=V3)) + #,y=log10(..count..))) + 
  geom_histogram(data=subset(mf2$DUP, mf2$DUP$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2$DUP, mf2$DUP$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Distance of CDS to Duplications") + #ylim(1,3) +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)



G4 <- ggplot(mf2$INV, aes(x=V3,y=log10(..count..))) + 
  geom_histogram(data=subset(mf2$INV, mf2$INV$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2$INV, mf2$INV$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Distance of CDS to Inversions") +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  #geom_density(aes(fill=factor(mf2$INV[,1])), alpha=0.15) +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)

G5 <- ggplot(mf2$TRA, aes(x=V3,y=log10(..count..))) + 
  geom_histogram(data=subset(mf2$TRA, mf2$TRA$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2$TRA, mf2$TRA$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Distance of CDS to Transpositions") +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=0, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)


#grid.arrange(G1,G2,G0, ncol=2)
grid.arrange(G1,G2,G3,G4,G5,G0, ncol=3)




dfm2<-read.table("signed_dist_ALL_Filt", header=F, sep = "\t")
mf2b <- split(dfm2, dfm2$V2)


G0b <- ggplot(dfm2, aes(x=V2, ..count..), fill=dfm2$V1) + scale_y_log10() + 
  geom_bar(aes(fill=dfm2$V1), position = "dodge") +
  ggtitle("Signed distance of CDS to Deletions") +
  geom_text(stat='count',aes(label=..count..),position = position_dodge(width = 0.5)) +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  theme_bw() + theme(legend.title=element_blank()) 


G1b <- ggplot(mf2b$DEL, aes(x=V3,y=log10(..count..))) + #scale_y_log10() + #scale_y_log10(limits = c(-0.1,NA)) + 
  geom_histogram(data=subset(mf2b$DEL, mf2b$DEL$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2b$DEL, mf2b$DEL$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Signed distance of CDS to Deletions") +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") + 
  theme_bw() + theme(legend.title=element_blank()) + #expand_limits(y=1) +
  scale_colour_manual(name="Selection-type", values=c(Positive="green", Non-positive="red")) +
  annotate("rect", xmin=1000, xmax=5000, ymin=-Inf, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)


G2b <- ggplot(mf2b$DUP, aes(x=V3)) + #,y=log10(..count..))) + 
  geom_histogram(data=subset(mf2b$DUP, mf2b$DUP$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2b$DUP, mf2b$DUP$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Signed distance of CDS to Duplications") + #ylim(1,3) +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=-Inf, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)

G3b <- ggplot(mf2b$INS, aes(x=V3,y=log10(..count..))) + 
  geom_histogram(data=subset(mf2b$INS, mf2b$INS$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2b$INS, mf2b$INS$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Signed distance of CDS to Insertions") +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=-Inf, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)

G4b <- ggplot(mf2b$INV, aes(x=V3,y=log10(..count..))) + 
  geom_histogram(data=subset(mf2b$INV, mf2b$INV$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2b$INV, mf2b$INV$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Signed distance of CDS to Inversions") +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  #geom_density(aes(fill=factor(mf2$INV[,1])), alpha=0.15) +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=-Inf, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)

G5b <- ggplot(mf2b$TRA, aes(x=V3,y=log10(..count..))) + 
  geom_histogram(data=subset(mf2b$TRA, mf2b$TRA$V1=='Positive'), aes(y=..count../sum(..count..)), fill="green", alpha=0.3) + 
  geom_histogram(data=subset(mf2b$TRA, mf2b$TRA$V1=='Non-positive'), aes(y=..count../sum(..count..)), fill="red", alpha=0.3) +
  ggtitle("Signed distance of CDS to Transpositions") +
  ylab("Frequency of distances") + xlab("Distance of SV to CDS feature") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate("rect", xmin=1000, xmax=5000, ymin=-Inf, ymax=Inf, fill="floralwhite",
           colour="brown", alpha=0.05)


grid.arrange(G1b,G2b,G3b,G4b,G5b,G0b, ncol=3)
