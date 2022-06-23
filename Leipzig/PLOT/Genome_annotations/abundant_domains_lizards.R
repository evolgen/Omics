library(ggplot2)
library(gridExtra)
#library(reshape2)

setwd("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/")

dat2 <- read.table(file = "/scr/bloodymary/rohit/Lacerta_viridis/Pfam/abundant_pfams.list",header = F,sep = "\t")

#dfm2 <- melt(dat2)[,c(1,2,4)]

svg(filename="/homes/biertank/rohit/Downloads/scripts/Genome_annotations/pfam_abundance_lizards.svg", 
    height = 8, width = 12, pointsize = 12)

#fam_names = c("ZF C2H2-type","Ankyrin repeats (multiple)","ZF double-domain","Immunoglobulin I-set",
#              "WD domain (WD40)","ZF C2H2 type","Ankyrin repeat","7 transmembrane receptor")

ggplot(dat2, aes(x=dat2$V1,y=log10(dat2$V3),fill=factor(dat2$V2))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  #xlab("pFAM domain names")+ 
  ylab("Abundance of pPFAM domains (log-scaled)") + 
#  labs(title="Most abundant protein domains in lizards") +
  scale_x_discrete(labels = c("PF00096.21" = "C2H2 type zinc finger","PF13637.1" = "Ankyrin repeats (multiple copies)",
              "PF13465.1" = "Zinc finger double domain", "PF07679.11" = "Immunoglobulin I-set",
              "PF00400.27" = "WD domain (WD40)", "PF13894.1" = "C2H2_4 type zinc finger ",
              "PF00023.25" = "Ankyrin repeat", "PF00001.16" = "7 transmembrane receptor")) + 
  # theme(legend.title=element_blank()) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(face = 'bold'),
                     axis.title.y = element_blank()) +
  scale_fill_manual(guide = guide_legend(title = "Species", reverse=TRUE),
                    values = c("orangered4","lightsalmon","blue2","green3")) + coord_flip()


dev.off()




