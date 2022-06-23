

# source("https://bioconductor.org/biocLite.R")
# biocLite("geneplotter")
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
library(gridExtra)
#library(xlsx)

setwd("/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts")
# Load directory and get mapping-files


#### Working version

directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_features")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_features",pattern="^new_features_")


### For Condition and Group
sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),5))
sampleCondition = c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,group=sampleGroup)

#ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group*condition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group+condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
colData(ddsHTSeq)$group<-factor(colData(ddsHTSeq)$group, levels=c('Bilineata','Viridis'))

ddsHTSeq <- ddsHTSeq[, ddsHTSeq$group %in% c('Bilineata','Viridis') ]
ddsHTSeq <- ddsHTSeq[, ddsHTSeq$condition %in% c('Brain','Gonads','Heart','Kidney','Liver') ]

dds<-DESeq(ddsHTSeq)
res<-results(dds)
summary(res)
res<-res[order(res$padj),]
head(res)

DESeq2::plotMA(dds,ylim=c(-11,11),main='DESeq2 with Groups and Conditions')
dev.copy(png,'~/Downloads/scripts/new_deseq2_MAplot_mod.png')
dev.off()

diff_genes <- res[which(res$pvalue < 0.05 & res$log2FoldChange > 1),]
write.table(diff_genes, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_differentially_expressed_mod.txt", sep="\t")

diff_genes_2 <- res[which(res$pvalue < 0.05 & res$log2FoldChange < 0),]
write.table(diff_genes_2, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_lowly_expressed_mod.txt", sep="\t")


ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)
summary(resLRT)
resLRT<-resLRT[order(res$padj),]
head(resLRT)

DESeq2::plotMA(ddsLRT,ylim=c(-11,11),main='DESeq2 with Groups and Conditions - LRT')
dev.copy(png,'~/Downloads/scripts/new_deseq2_MAplot_mod_LRT.png')
dev.off()

diff_genes_LRT <- resLRT[which(resLRT$pvalue < 0.05 & resLRT$log2FoldChange > 1),]
write.table(diff_genes_LRT, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_differentially_expressed_mod_LRT.txt", sep="\t")

diff_genes_2_LRT <- resLRT[which(resLRT$pvalue < 0.05 & resLRT$log2FoldChange < 0),]
write.table(diff_genes_2_LRT, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_lowly_expressed_mod_LRT.txt", sep="\t")


### For Groups only
directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_features")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_features",pattern="^new_features_")

sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),5))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, group=sampleGroup)
ddsHTSeq_group<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group)
colData(ddsHTSeq_group)$group<-factor(colData(ddsHTSeq_group)$group, levels=c('Bilineata','Viridis'))

ddsHTSeq_group <- ddsHTSeq_group[, ddsHTSeq_group$group %in% c('Bilineata','Viridis') ]

dds_group<-DESeq(ddsHTSeq_group)
res_group<-results(dds_group)
res_group<-res_group[order(res_group$padj),]
head(res_group)

DESeq2::plotMA(dds_group,ylim=c(-5,5),main='DESeq2 with Groups only')
dev.copy(png,'~/Downloads/scripts/new_deseq2_MAplot_group.png')
dev.off()

diff_genes_group <- res_group[which(res_group$pvalue < 0.05 & res_group$log2FoldChange > 1),]
write.table(diff_genes_group, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_differentially_expressed_group.txt", sep="\t")

diff_genes_group_2 <- res_group[which(res_group$pvalue < 0.05 & res_group$log2FoldChange < 0),]
write.table(diff_genes_group_2, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_lowly_expressed_group.txt", sep="\t")



### For Conditions only
directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_features")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_features",pattern="^new_features_")

sampleCondition = c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

ddsHTSeq_cond<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq_cond)$condition<-factor(colData(ddsHTSeq_cond)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
ddsHTSeq_cond <- ddsHTSeq_cond[, ddsHTSeq_cond$condition %in% c('Brain','Gonads','Heart','Kidney','Liver') ]

dds_cond<-DESeq(ddsHTSeq_cond)
res_cond<-results(dds_cond)
res_cond<-res_cond[order(res_cond$padj),]
head(res_cond)

dds_cond_LRT <- DESeq(ddsHTSeq_cond, test="LRT", reduced= ~ 1)
res_cond_LRT <- results(dds_cond_LRT)

DESeq2::plotMA(dds_cond,ylim=c(-11,11),main='DESeq2 with Conditions only')
dev.copy(png,'~/Downloads/scripts/new_deseq2_MAplot_cond.png')
dev.off()

diff_genes_cond <- res_cond[which(res_cond$pvalue < 0.05 & res_cond$log2FoldChange > 1),]
write.table(diff_genes_cond, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_differentially_expressed_cond.txt", sep="\t")

diff_genes_cond_2 <- res_cond[which(res_cond$pvalue < 0.05 & res_cond$log2FoldChange < 0),]
write.table(diff_genes_cond_2, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_lowly_expressed_cond.txt", sep="\t")




dds_cond<-DESeq(ddsHTSeq_cond)
res.B.G <- results(dds_cond, contrast=c("condition","Brain","Gonads"))
res.B.G <-res.B.G [order(res.B.G $padj),]
head(res.B.G)
res.B.H <- results(dds_cond, contrast=c("condition","Brain","Heart"))
res.B.H<-res.B.H[order(res.B.H$padj),]
res.B.L <- results(dds_cond, contrast=c("condition","Brain","Liver"))
res.B.L<-res.B.L[order(res.B.L$padj),]
res.B.K <- results(dds_cond, contrast=c("condition","Brain","Kidney"))
res.B.K<-res.B.K[order(res.B.K$padj),]

res.G.H <- results(dds_cond, contrast=c("condition","Gonads","Heart"))
res.G.H<-res.G.H[order(res.G.H$padj),]
res.G.K <- results(dds_cond, contrast=c("condition","Gonads","Kidney"))
res.G.K<-res.G.K[order(res.G.K$padj),]
res.G.L <- results(dds_cond, contrast=c("condition","Gonads","Liver"))
res.G.L<-res.G.L[order(res.G.L$padj),]

res.H.K <- results(dds_cond, contrast=c("condition","Heart","Kidney"))
res.H.K<-res.H.K[order(res.H.K$padj),]
res.H.L <- results(dds_cond, contrast=c("condition","Heart","Liver"))
res.H.L<-res.H.L[order(res.H.L$padj),]

res.K.L <- results(dds_cond, contrast=c("condition","Kidney","Liver"))
res.K.L<-res.K.L[order(res.K.L$padj),]


par(mfrow=c(5,2))
DESeq2::plotMA(res.B.G,ylim=c(-11,11),main='DESeq2 with Conditions - Brain vs. Gonads')
DESeq2::plotMA(res.B.H,ylim=c(-11,11),main='DESeq2 with Conditions - Brain vs. Heart')
DESeq2::plotMA(res.B.L,ylim=c(-11,11),main='DESeq2 with Conditions - Brain vs. Liver')
DESeq2::plotMA(res.B.K,ylim=c(-11,11),main='DESeq2 with Conditions - Brain vs. Kidney')
DESeq2::plotMA(res.G.H,ylim=c(-11,11),main='DESeq2 with Conditions - Gonads vs. Heart')
DESeq2::plotMA(res.G.L,ylim=c(-11,11),main='DESeq2 with Conditions - Gonads vs. Liver')
DESeq2::plotMA(res.G.K,ylim=c(-11,11),main='DESeq2 with Conditions - Gonads vs. Kidney')
DESeq2::plotMA(res.H.L,ylim=c(-11,11),main='DESeq2 with Conditions - Heart vs. Liver')
DESeq2::plotMA(res.H.K,ylim=c(-11,11),main='DESeq2 with Conditions - Heart vs. Kidney')
DESeq2::plotMA(res.K.L,ylim=c(-11,11),main='DESeq2 with Conditions - Kidney vs. Liver')
dev.copy(png,'~/Downloads/scripts/new_deseq2_MAplot_cond_pairwise.png')
dev.off()


