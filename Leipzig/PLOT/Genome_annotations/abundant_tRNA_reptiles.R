library(ggplot2)
library(gridExtra)
#library(reshape2)

setwd("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/")

dat2 <- read.table(file = "miR",header = T,sep = "\t")

dfm2 <- melt(dat2)[,c(1,2,4)]


bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/tRNA_abundance_infernal.tiff", height = 8, width = 12, units = 'in', res=1600)
 

# ggplot(dfm2,aes(x=tRNA,y=value,fill=factor(Species))) +
#   geom_bar(stat="identity",position="dodge") + 
#   xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + labs(title="tRNAs abundant in Lacertids") +
#   theme_bw() + theme(legend.title=element_blank()) # + scale_fill_manual(values = c("blue2","green3")) + 
#   #geom_text(aes(label=value), vjust=-0.5, position=position_dodge(width=1), size=4)


trnaplot1 <- ggplot(subset(dfm2,tRNA %in% c("Asn_GTT" , "Glu_CTC", "Lys_TTT", "Pseudo_ATT", "Pseudo_CAT", "Pseudo_GTT", "Pseudo_CCT", "Pseudo_CTC", "Pseudo_TTT")), 
       aes(x=tRNA,y=log(value),fill=factor(Species))) +
  geom_bar(stat="identity",position="dodge", colour="black") + 
  xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + 
  labs(title="Highly abundant tRNAs in Lacertids - I") +
  theme_bw() + theme(legend.title=element_blank()) + 
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","blue2","green3")) # + 
#geom_text(aes(label=value), vjust=-0.5, position=position_dodge(width=1), size=4)

trnaplot2 <- ggplot(subset(dfm2,!(tRNA %in% c("Asn_GTT" , "Glu_CTC", "Lys_TTT", "Pseudo_ATT", "Pseudo_CAT", "Pseudo_GTT", "Pseudo_CCT", "Pseudo_CTC", "Pseudo_TTT"))), 
       aes(x=tRNA,y=log(value),fill=factor(Species))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + 
  labs(title="Highly abundant tRNAs in Lacertids - II") +
  theme_bw() + theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","blue2","green3")) 


grid.arrange(trnaplot1, trnaplot2, ncol=1)

dev.off()




