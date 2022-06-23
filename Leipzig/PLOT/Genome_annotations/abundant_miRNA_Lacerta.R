library(ggplot2)
library(gridExtra)
#library(reshape2)

setwd("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/")

dat2 <- read.table(file = "miRNA_infernalonly.txt",header = F,sep = "\t")

#dfm2 <- melt(dat2)[,c(1,2,4)]


svg(filename="/homes/biertank/rohit/Downloads/scripts/Genome_annotations/miRNA_abundance_infernal.svg", height = 8, width = 12, pointsize = 12)


ggplot(dat2, aes(x=dat2$V3,y=log10(dat2$V2),fill=factor(dat2$V1))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  xlab("miRNA family")+ylab("Abundance of miRNA (log-scaled)") + 
  labs(title="Abundant miRNAs detected by Infernal without blast homologs") +
  theme_bw() + theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("blue2","green3")) 


dev.off()




