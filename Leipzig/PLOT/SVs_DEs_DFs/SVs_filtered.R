setwd("C:/Users/rohit/Desktop/GenomePaper/SVs/")

library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(reshape2)

dfm2<-read.table("SVs_filtered.txt", header=TRUE, sep = "\t")
#dfm2<-read.table("filterd_SVs_ranges.txt", header=TRUE, sep = "\t")


df.1<-melt(dfm2, id.vars=c(1))

cols <- colorRampPalette(brewer.pal(5, "Set1"))
ngroups <- length(unique(df.1$variable))

#orders <- c("Complex"="Complex","DEL"="Deletion","IDP"="Break-point",
#            "INS"="Insertion","INV"="Inversion","INVDUP"="INV+DUP","TRA"="Translocation")

svg("/Users/rohit/Desktop/GenomePaper/SVs/SVs_filtered.svg", 
              height = 8, width = 12, pointsize=12)
       
ggplot(df.1,aes(x=Type,y=log10(value),fill=variable)) +
  geom_histogram(stat="identity",position="dodge") + theme_bw() +
  ggtitle("Final list of SVs detected in Lacertids") +
  scale_fill_discrete(name ="Type of variable") +
  ylab("Size frequency of SVs (log-scaled)") + xlab("Type of SV") +
  theme_bw() + scale_fill_manual(values = cols(ngroups)) +
  theme_bw() + theme(legend.title=element_blank())

dev.off()
