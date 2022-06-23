
#source("https://bioconductor.org/biocLite.R")

library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(RcppArmadillo)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(stringr)
library(genefilter)
library(biomaRt)
library(dplyr)
library(fdrtool)
library(gridExtra)
library(data.table)
library(tidyr)
library(vsn)
library(ggpubr)


setwd("/scr/k70san/rohit/Lacertidae/splicing/pooled/ALL/HTseq/TRIAD/")

#### Working version

directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/pooled/ALL/HTseq/TRIAD/Counts")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/pooled/ALL/HTseq/TRIAD/Counts",pattern="^new_features_")
dds<-read.table("/scr/k70san/rohit/Lacertidae/splicing/pooled/ALL/HTseq/TRIAD/Counts/SampleTable.txt.extended.edit", 
                header=T, sep = "\t")

dds<-DESeqDataSetFromHTSeqCount(sampleTable=dds, directory=directory, design=~group*condition)

dds <- dds[ rowSums(counts(dds)) > 100, ]
rld <- rlog(dds, blind=FALSE)      #For visualising
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
ntd <- normTransform(dds)
normalised_counts <- assay(ntd)
write.table(normalised_counts, file='Normalized_readcounts.txt')


p1 <- plotPCA(vsd, intgroup=c("group"))
p2 <- plotPCA(vsd, intgroup=c("condition"))
p3 <- plotPCA(vsd, intgroup=c("sex"))
p4 <- plotPCA(vsd, intgroup=c("group","condition","sex"))
ggarrange(p1, p2, p3,p4 + rremove("x.text"), labels = c("group", "condition", "group+condition"), ncol = 2, nrow = 2)


plotrld <- plotPCA(vsd, intgroup=c("group")) + theme_bw() + labs(col="Tissue")
plotvsd <- plotPCA(rld, intgroup=c("group")) + theme_bw() + labs(col="Tissue")
svg(filename="/homes/biertank/rohit/Downloads/scripts/differential_expression/tissue_specific_view.svg", 
    height = 8, width = 12, pointsize = 12)  
grid.arrange(grobs = list(plotrld, plotvsd), ncol=2)
dev.off()


plotrld <- plotPCA(vsd, intgroup=c("group","condition","sex")) + theme_bw() + labs(col="Tissue+Species") +  
  theme(legend.title=element_blank(), axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14,face="bold"))
plotvsd <- plotPCA(rld, intgroup=c("group","condition","sex")) + theme_bw() + labs(col="Tissue+Species") +
  theme(legend.title=element_blank(), axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=14,face="bold"))
svg(filename="/homes/biertank/rohit/Downloads/scripts/differential_expression/tissue_specific_view_detailed.svg", 
    height = 8, width = 12, pointsize = 12)  
grid.arrange(grobs = list(plotrld, plotvsd), ncol=2)
dev.off()


#dds <- DESeq(dds)
#resultsNames(dds)

