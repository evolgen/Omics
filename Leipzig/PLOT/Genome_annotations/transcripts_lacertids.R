library(ggplot2)
library(gridExtra)

setwd("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/")

dat2 <- read.table(file = "/homes/biertank/rohit/Downloads/scripts/Genome_annotations/transcripts_lacertids.txt",header = F,sep = "\t")


svg(filename="/homes/biertank/rohit/Downloads/scripts/Genome_annotations/transcripts_paper.svg", 
    height = 8, width = 12, pointsize = 12)


ggplot(dat2, aes(x=dat2$V1,y=(dat2$V2),fill=factor(dat2$V3))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  xlab("Tissue source of the transcripts") + ylab("Number of transcripts (log-scaled)") + 
  #labs(title="Segmental duplications in Lacertids") + 
  coord_cartesian(ylim=c(0,60000)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(face = 'bold'),
                     axis.title.y = element_text('bold')) +
  scale_fill_manual(values = c("blue2","green3")) + theme(legend.title=element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12, face = "bold"),
        legend.text=element_text(size=13, face = "bold"), legend.key.size = unit(1.25,"line"))



dev.off()
