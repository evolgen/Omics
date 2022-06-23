

source("https://bioconductor.org/biocLite.R")
biocLite("geneplotter")
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(RcppArmadillo)
#library(DESeq)
 library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
#library(xlsx)

setwd("/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts")
# Load directory and get mapping-files


#### Working version

directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts",pattern="^features_")

sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),5))
sampleCondition = c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,group=sampleGroup)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group*condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
colData(ddsHTSeq)$group<-factor(colData(ddsHTSeq)$group, levels=c('Bilineata','Viridis'))


ddsHTSeq <- ddsHTSeq[, ddsHTSeq$group %in% c('Bilineata','Viridis') ]
ddsHTSeq <- ddsHTSeq[, ddsHTSeq$condition %in% c('Brain','Gonads','Heart','Kidney','Liver') ]

# ddsHTSeq$condition <- droplevels(ddsHTSeq$condition)

dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)

DESeq2::plotMA(dds,ylim=c(-2,2),main='DESeq2')
dev.copy(png,'~/Downloads/scripts/deseq2_MAplot.png')
dev.off()

diff_genes <- res[which(res$pvalue < 0.05 & res$log2FoldChange > 1),]
write.table(diff_genes, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/differentially_expressed.txt", sep="\t")

diff_genes_2 <- res[which(res$pvalue < 0.05 & res$log2FoldChange < 0),]
write.table(diff_genes_2, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/lowly_expressed.txt", sep="\t")


#dds_norepl = dds[ , c( "features_brain_cross.txt", "features_brain_self.txt" ) ]
#dds_norepl = estimateDispersions( dds_norepl, method="blind", sharingMode="fit-only" )
#results_2 = nbinomTest(dds_norepl, "Bilineata", "Viridis")



