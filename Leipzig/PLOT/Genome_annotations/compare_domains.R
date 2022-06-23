library(ggplot2)
library(gridExtra)

setwd("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/")

dat1 <- read.table(file = "/scr/bloodymary/rohit/Lacerta_viridis/Interpro/domain_pfam_lacertids_top6.txt",header = F,sep = "\t")

plot2 <- ggplot(dat1, aes(x=dat1$V2,y=dat1$V3,fill=factor(dat1$V4))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  ylab("Number of pPFAM domains (log-scaled)") + 
  labs(title="(B)                                            ") +
  scale_x_discrete(labels = c("PF00069" = "Protein kinase domain", "PF00076" = "RNA recognition motif (RNP domain)",
                              "PF00096" = "Zinc finger, C2H2 type", "PF00400" = "WD domain, G-beta repeat",
                              "PF07679" = "Immunoglobulin I-set domain", "PF12796" = "Ankyrin repeats (3 copies)")) +
  theme_bw() + theme(axis.text=element_text(size=14,face="bold"),
                     axis.title=element_text(size=14,face="bold"),
                     axis.title.y = element_blank(),legend.position="none") +
  scale_fill_manual(guide = guide_legend(title = "Species", reverse=TRUE),
                    values = c("blue2","green3")) + coord_flip()
  

dat2 <- read.table(file = "/scr/bloodymary/rohit/Lacerta_viridis/Pfam/abundant_pfams.list",header = F,sep = "\t")

plot1 <- ggplot(dat2, aes(x=dat2$V1,y=log10(dat2$V3),fill=factor(dat2$V2))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  ylab("Number of gene family categories (log-scaled)") + 
  labs(title="(A)                                            ") +
  scale_x_discrete(labels = c("PF00096.21" = "C2H2 type zinc finger","PF13637.1" = "Ankyrin repeats (multiple copies)",
                              "PF13465.1" = "Zinc finger double domain", "PF07679.11" = "Immunoglobulin I-set",
                              "PF00400.27" = "WD domain (WD40)", "PF13894.1" = "C2H2_4 type zinc finger ",
                              "PF00023.25" = "Ankyrin repeat", "PF00001.16" = "7 transmembrane receptor")) + 
  theme_bw() + theme(axis.text=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"), 
                     axis.title.y = element_blank(), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14,face="bold.italic")) +
  scale_fill_manual(guide = guide_legend(title = "Species", reverse=TRUE),
                    values = c("orangered4","lightsalmon","blue2","green3")) + coord_flip()


svg(filename="/homes/biertank/rohit/Downloads/scripts/Genome_annotations/protein_domains_lizards_lacertids.svg", 
    height = 8, width = 12, pointsize = 12)

grid.arrange(grobs = list(plot1, plot2), ncol=1)

dev.off()
