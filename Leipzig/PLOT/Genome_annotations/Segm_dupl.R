library(ggplot2)
library(gridExtra)

setwd("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/")

dat2 <- read.table(file = "/homes/biertank/rohit/Downloads/scripts/Genome_annotations/Segm_dupl_lacertids.txt",header = F,sep = "\t")


svg(filename="/homes/biertank/rohit/Downloads/scripts/Genome_annotations/segmental_duplications_paper.svg", 
    height = 8, width = 12, pointsize = 12)


ggplot(dat2, aes(x=dat2$V1,y=log10(dat2$V2),fill=factor(dat2$V3))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  xlab("Minimum size of segmental duplication") + ylab("Frequency of segmental duplications (log-scaled)") + 
  #labs(title="Segmental duplications in Lacertids") + 
  scale_x_continuous(breaks = c(1000,5000,10000)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(face = 'bold'),
                     axis.title.y = element_text('bold')) +
  scale_fill_manual(values = c("blue2","green3")) + theme(legend.title=element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12, face = "bold"),
        legend.text=element_text(size=13, face = "bold"), legend.key.size = unit(1.25,"line"))



dev.off()
