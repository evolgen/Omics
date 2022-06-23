

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
library(pheatmap)
library(PoiClaClu)
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

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group+condition+group:condition)
#ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group*condition)
#ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group+condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
colData(ddsHTSeq)$group<-factor(colData(ddsHTSeq)$group, levels=c('Bilineata','Viridis'))

ddsHTSeq <- ddsHTSeq[, ddsHTSeq$group %in% c('Bilineata','Viridis') ]
ddsHTSeq <- ddsHTSeq[, ddsHTSeq$condition %in% c('Brain','Gonads','Heart','Kidney','Liver') ]

dds<-DESeq(ddsHTSeq)
res<-results(dds)
summary(res)
res<-res[order(res$padj),]
head(res)
resultsNames(dds)

countdata<-assay(ddsHTSeq)
coldata <- colData(ddsHTSeq)

nrow(ddsHTSeq)
#33824
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]
#33797
rld <- rlog(ddsHTSeq, blind=TRUE)
head(assay(rld), 3)

png(file="~/Downloads/scripts/differential_expression/new_deseq2_scatterplot_counts.png",width=1200,height=1100)
par( mfrow = c( 1, 2 ) )
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
plot(log2(counts(ddsHTSeq, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.3)
plot(assay(rld)[,1:2], pch=16, cex=0.3)
# plot(assay(rld)[,3:4], pch=16, cex=0.3)
# plot(assay(rld)[,5:6], pch=16, cex=0.3)
# plot(assay(rld)[,7:8], pch=16, cex=0.3)
# plot(assay(rld)[,9:10], pch=16, cex=0.3)
dev.off()

sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$group, rld$condition, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(file="~/Downloads/scripts/differential_expression/new_deseq2_heatplot_distmat.png",width=1200,height=1100)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off()

poisd <- PoissonDistance(t(counts(ddsHTSeq)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$group, rld$condition, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
png(file="~/Downloads/scripts/differential_expression/new_deseq2_heatplot_poisdist.png",width=1200,height=1100)
pheatmap(samplePoisDistMatrix, clustering_distance_rows=poisd$dd, clustering_distance_cols=poisd$dd, col=colors)
dev.off()

png(file="~/Downloads/scripts/differential_expression/new_deseq2_PCA.png",width=1200,height=1100)
plotPCA(rld, intgroup = c("group","condition"))
dev.off()

data <- plotPCA(rld, intgroup = c("group","condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

# par( mfrow = c( 1, 2 ) )
# p1 <- pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
# p2 <- pheatmap(samplePoisDistMatrix, clustering_distance_rows=poisd$dd, clustering_distance_cols=poisd$dd, col=colors)


### MA-plot
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
resultsNames(dds_group)


DESeq2::plotMA(dds_group,ylim=c(-5,5),main='DESeq2 with Groups only')
dev.copy(png,'~/Downloads/scripts/new_deseq2_MAplot_group.png')
dev.off()

diff_genes_group <- res_group[which(res_group$pvalue < 0.05 & res_group$log2FoldChange > 1),]
write.table(diff_genes_group, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_differentially_expressed_group.txt", sep="\t")

diff_genes_group_2 <- res_group[which(res_group$pvalue < 0.05 & res_group$log2FoldChange < 0),]
write.table(diff_genes_group_2, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_lowly_expressed_group.txt", sep="\t")



### For conditions only
ddsHTSeq_condition<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq_condition)$group<-factor(colData(ddsHTSeq_condition)$group, levels=c('Bilineata','Viridis'))

ddsHTSeq_condition <- ddsHTSeq_condition[, ddsHTSeq_condition$group %in% c('Bilineata','Viridis') ]

dds_condition<-DESeq(ddsHTSeq_condition)
res_condition<-results(dds_condition)
res_condition<-res_condition[order(res_condition$padj),]
head(res_condition)

DESeq2::plotMA(dds_condition,ylim=c(-5,5),main='DESeq2 with Conditions only')
dev.copy(png,'~/Downloads/scripts/new_deseq2_MAplot_condition.png')
dev.off()

diff_genes_condition <- res_condition[which(res_condition$pvalue < 0.05 & res_group$log2FoldChange > 1),]
write.table(diff_genes_condition, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_differentially_expressed_condition.txt", sep="\t")

diff_genes_condition_2 <- res_condition[which(res_condition$pvalue < 0.05 & res_group$log2FoldChange < 0),]
write.table(diff_genes_condition_2, "/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_lowly_expressed_condition.txt", sep="\t")


results(dds_condition, contrast = list(c("Brain"), c("Gonads","Heart","Kidney","Liver")), listValues=c(1,-1/4))



### For conditions affecting groups
ddsHTSeq_condition_G<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition+group)
colData(ddsHTSeq_condition_G)$group<-factor(colData(ddsHTSeq_condition_G)$group, levels=c('Bilineata','Viridis'))

ddsHTSeq_condition_G <- ddsHTSeq_condition_G[, ddsHTSeq_condition_G$group %in% c('Bilineata','Viridis') ]

dds_condition_G<-DESeq(ddsHTSeq_condition_G)
res_condition_G<-results(dds_condition_G)
res_condition_G<-res_condition_G[order(res_condition_G$padj),]
head(res_condition_G)

resultsNames(dds_condition_G)

results(dds_condition_G, contrast = list(c("Brain"), c("Gonads","Heart","Kidney","Liver")), listValues=c(1,-1/4))



### Combination
directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/pooled/namesort/counts/new_features")
sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),5))
sampleCondition = c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,group=sampleGroup)

dds <- sampleTable
dds$group2 <- factor(paste0(dds$condition, dds$group))

dds<-DESeqDataSetFromHTSeqCount(sampleTable=dds, 
                                     directory=directory, design=~group2)

#design(dds) <- ~ group2
dds <- DESeq(dds)
resultsNames(dds)
# e.g. for condition KO cell type 2 vs cell type 1
results(dds, contrast=c("group2","GonadsBilineata","GonadsViridis")) 






