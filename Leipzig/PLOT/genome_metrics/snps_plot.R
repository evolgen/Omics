 rm(list=ls())
library(psych)
library(binom)
z<-
read.table("short_introg_indels_edited_plot")
x<-
read.table("short_introg_snps")
u<-rbind(z$V1/5430,x$V1/146816) ### number of del and ins, normliasing
colnames(u)<-z$V2
rownames(u)<-c("introgressed indels","introgressed SNVs")
png("new_vep_intro_indels_snps.png",
width = 10, height = 6, units = 'in', res = 400)
barplot(u,beside=T,legend.text=T,col=c("green","darkorange1"),xlab="",
ylab="Fraction", main="Annotation of Introgressed Indel and
SNP's",cex.names=0.5,las=2,axisnames = FALSE, ylim = c(0,0.65))
text(seq(2.5,ncol(u)*3,3), -0.01, srt = 15, adj= 1, xpd = TRUE,labels =
colnames(u), cex=0.7)

conu1<-binom.confint(z$V1,5430, conf.level = 0.95, methods = "exact") # for
introgressed and shared indels
conu2<-binom.confint(x$V1,146816,conf.level = 0.95, methods = "exact")
aa<-seq(1.5,NROW(conu1)*3,3)
bb<-seq(2.5,NROW(conu2)*3,3)
segments(aa,conu1$lower,x1=aa,y1=conu1$upper,lwd=1)
segments(aa-0.1,conu1$lower,aa+0.1,y1=conu1$lower,lwd=1)
segments(aa-0.1,conu1$upper,aa+0.1,y1=conu1$upper,lwd=1)
segments(bb,conu2$lower,x1=bb,y1=conu2$upper,lwd=1)
segments(bb-0.1,conu2$lower,bb+0.1,y1=conu2$lower,lwd=1)
segments(bb-0.1,conu2$upper,bb+0.1,y1=conu2$upper,lwd=1)
dev.off()


