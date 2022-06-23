library(ggplot2)
library(gridExtra)
library(reshape)
require(grid)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/Noncoding/tRNA/")

dat2 <- read.table(file = "abundant_lacertids_allcount.txt",header = T,sep = "\t")

dfm2 <- melt(dat2)[,c(1,2,4)]

dfm2 <- subset(dfm2,Species %in% c("lvi" , "lbi", "hsa", "gja", "acr", "ams", "cmy", "pbi"))

dfm2$Species <- factor(dfm2$Species, levels = c("cmy", "ams", "pbi", "acr", "gja", "lvi", "lbi", "hsa"))
dfm2 <- dfm2[order(dfm2$Species), ]

spec_order <- c("cmy", "ams", "pbi", "acr", "gja", "lvi", "lbi", "hsa")
dfm2$Species <- ordered(dfm2$Species, spec_order)
dfm2 <- dfm2[with(dfm2, order(tRNA, Species)),]


bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/tRNA_abundance_lacertids_3.tiff", height = 12, width = 21, units = 'in', res=1600)


types <- unique(dfm2$tRNA)
type_1 <- types[1:7]
type_2 <- types[8:14]
type_3 <- types[15:21]
type_4 <- types[22:28]
type_5 <- types[29:36]
type_6 <- types[37:44]

# ggplot(dfm2,aes(x=tRNA,y=value,fill=factor(Species))) +
#   geom_bar(stat="identity",position="dodge") + 
#   xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log)") + labs(title="tRNAs abundant in Lacertids") +
#   theme_bw() + theme(legend.title=element_blank()) # + scale_fill_manual(values = c("blue2","green3")) + 
#   #geom_text(aes(label=value), vjust=-0.5, position=position_dodge(width=1), size=4)


trnaplot1 <- ggplot(subset(dfm2,tRNA %in% type_1), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity",position="dodge", colour="black") + 
  #xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + 
  #labs(title="Highly abundant tRNAs in Lacertids - I") +
  theme_bw() + theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"))
#+ geom_text(aes(label=value), vjust=-0.5, position=position_dodge(width=1), size=4)

trnaplot2 <- ggplot(subset(dfm2,tRNA %in% type_2), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  #xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + 
  #labs(title="Highly abundant tRNAs in Lacertids - II") +
  theme_bw() + theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"))
#+ scale_fill_manual(values = c("firebrick3","olxivedrab","tan4","blue2","green3")) 

trnaplot3 <- ggplot(subset(dfm2,tRNA %in% type_3), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  #xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + 
  #labs(title="Highly abundant tRNAs in Lacertids - III") +
  theme_bw() + theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"))


trnaplot4 <- ggplot(subset(dfm2,tRNA %in% type_4), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  #xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + 
  #labs(title="Highly abundant tRNAs in Lacertids - IV") +
  theme_bw() + theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"))

trnaplot5 <- ggplot(subset(dfm2,tRNA %in% type_5), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  xlab("Type of the tRNA with anti-codon")+ #ylab("Abundance of tRNA (log-scaled)") + 
  #labs(title="Highly abundant tRNAs in Lacertids - V") +
  theme_bw() + theme(legend.title=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"))

trnaplot6 <- ggplot(subset(dfm2,tRNA %in% type_6), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  xlab("Type of the tRNA with anti-codon")+ #ylab("Abundance of tRNA (log-scaled)") + 
  #labs(title="Highly abundant tRNAs in Lacertids - VI") +
  theme_bw() + theme(legend.title=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"))


grid.arrange(trnaplot1, trnaplot2, trnaplot3, trnaplot4, trnaplot5, trnaplot6, ncol=2,
             top=textGrob("Highly abundant tRNAs in Lacertids (tRNA-SCAN)", gp=gpar(fontsize=16,font=8)))

dev.off()


