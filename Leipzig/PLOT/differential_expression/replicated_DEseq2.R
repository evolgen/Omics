

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

setwd("/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts")
# Load directory and get mapping-files


#### Working version

directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts",pattern="^new_features_")


### For Condition and Group
sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),15))
sampleCondition = c(rep.int(c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2)),3))
sampleReplicate <- c(rep("1",10), rep("2",10),rep("3",10))
#sampleReplicate <- factor(rep.int(c("1","2","3"),5))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,
                        group=sampleGroup, replicate=sampleReplicate)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design =~ group*condition)
#ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design =~ group+condition+group:condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
colData(ddsHTSeq)$group<-factor(colData(ddsHTSeq)$group, levels=c('Bilineata','Viridis'))
#colData(ddsHTSeq)$replicate<-factor(colData(ddsHTSeq)$replicate, levels=c('1','2','3'))

ddsHTSeq <- ddsHTSeq[, ddsHTSeq$group %in% c('Bilineata','Viridis') ]
ddsHTSeq <- ddsHTSeq[, ddsHTSeq$condition %in% c('Brain','Gonads','Heart','Kidney','Liver') ]
#ddsHTSeq <- ddsHTSeq[, ddsHTSeq$replicate %in% c('1','2','3') ]
ddsHTSeq$group <- relevel(ddsHTSeq$group, "Viridis")

nrow(ddsHTSeq)
#33824
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]
#33780
rld <- rlog(ddsHTSeq, blind=TRUE)
head(assay(rld), 3)


## To detect outliers
plotPCA(ddsHTSeq, intgroup = c("group","condition"))
data <- plotPCA(rld, intgroup = c("group","condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
png(file="~/Downloads/scripts/differential_expression/reps_deseq2_PCA.png",width=1200,height=1100)
ggplot(data, aes(PC1, PC2, color=group.1, shape=condition)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()
dev.off()

outliers <- as.character(subset(colnames(ddsHTSeq), data$PC1 > 0))
outliers

dds <- estimateSizeFactors(ddsHTSeq)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
#idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

GeneCounts <- counts(ddsHTSeq)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)
#15425


## To check if normalisation worked with multidensity and multiecdf plots, pairwise MDplots
multidensity( counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_multidensity_mod.png')
dev.off()
multiecdf( counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1000))
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_multiecdf_mod.png')
dev.off()

#png(file="~/Downloads/scripts/differential_expression/reps_deseq2_MDplot_pairwise.png",width=1200,height=1100)
pdf("~/Downloads/scripts/differential_expression/reps_deseq2_MDplot_pairwise.pdf")
MA.idx = t(combn(1:30, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,], 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(ddsHTSeq)[MA.idx[i,1]], " vs ",
                       colnames(ddsHTSeq)[MA.idx[i,2]] ), ylim = c(-3,3))
}
dev.off()


png(file="~/Downloads/scripts/differential_expression/reps_deseq2_dispest.png",width=1200,height=1100)
plotDispEsts(estimateDispersions(estimateSizeFactors(ddsHTSeq)))
dev.off()

#resultsNames(dds)
#res<-results(dds,name="condition_Gonads_vs_Brain")
res<-results(dds,contrast=c("group","Bilineata","Viridis"))
res_Group_Gonads<-results(dds,contrast=list(c("group_Bilineata_vs_Viridis","groupBilineata.conditionGonads")))
summary(res_Group_Gonads)
# out of 32610 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 6419, 20% 
# LFC < 0 (down)   : 4945, 15% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 1)
res2G<-results(dds,list(c("groupBilineata.conditionGonads")))
summary(res2G)
res2G<-res2G[order(res$padj),]
# out of 32610 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 5680, 17% 
# LFC < 0 (down)   : 4536, 14% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 633, 1.9% 
# (mean count < 2)

#res<-results(dds, list( c("Intercept") ))
summary(res)
### Interaction * with zero-count removal (>=5)
# out of 32610 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 5243, 16% 
# LFC < 0 (down)   : 5070, 16% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 633, 1.9% 
# (mean count < 2)
### Interaction * with zero-count removal (>=15)
# out of 32610 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 5174, 16% 
# LFC < 0 (down)   : 5562, 17% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 633, 1.9% 
### Interaction * with zero-count removal (>=100)
# out of 11732 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 3494, 30% 
# LFC < 0 (down)   : 3241, 28% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 12)
### Interaction * without zero-count removal
# out of 33800 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 5104, 15% 
# LFC < 0 (down)   : 5495, 16% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 1311, 3.9% 
### Interaction + without zero-count removal
# out of 33800 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 8644, 26% 
# LFC < 0 (down)   : 11016, 33% 
# outliers [1]     : 58, 0.17% 
# low counts [2]   : 0, 0%
res<-res[order(res$padj),]
head(res)


plotMA(dds,ylim=c(-11,11),main='DESeq2 with Groups and Conditions')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_mod.png')
dev.off()

resLFC1 <- results(dds, lfcThreshold=1)
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="darkblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="darkblue")
})


## Check for filtered-ones, pvalue and correct them accordingly
ddsHTSeq_2 <-  nbinomWaldTest(estimateDispersions(estimateSizeFactors(ddsHTSeq)))
res_2 <- results(ddsHTSeq_2, pAdjustMethod = "BH",contrast=c("group","Bilineata","Viridis"))
res_2G <- results(ddsHTSeq_2, pAdjustMethod = "BH",contrast=list(c("groupBilineata.conditionGonads")))
res_Group_Gonads <- results(ddsHTSeq_2, pAdjustMethod = "BH",
                            contrast=list(c("group_Bilineata_vs_Viridis","groupBilineata.conditionGonads")))
table(res_2$padj < 0.1)
# FALSE  TRUE 
# 22980 10145 


png(file="~/Downloads/scripts/differential_expression/reps_deseq2_filterrej.png",width=1200,height=1100)
par(mfrow=c(1,2))
plot(metadata(res_2)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")
hist(res_2$pvalue, col = "lavender", main = "WT vs Deletion", xlab = "p-values")
dev.off()

res_2 <- res_2[ !is.na(res_2$padj), ]
res_2 <- res_2[ !is.na(res_2$pvalue), ]
res_2 <- res_2[, -which(names(res_2) == "padj")]
png(file="~/Downloads/scripts/differential_expression/reps_deseq2_pvaluecheck.png",width=1200,height=1100)
FDR.res_2 <- fdrtool(res_2$stat, statistic= "normal", plot = T)
dev.off()
FDR.res_2$param[1, "sd"]
# sd 
# 1.852826

png(file="~/Downloads/scripts/differential_expression/reps_deseq2_pvaluecorrected.png",width=1200,height=1100)
res_2[,"padj"]  <- p.adjust(FDR.res_2$pval, method = "BH")
hist(FDR.res_2$pval, col = "royalblue4", 
     main = "Viridis vs Bilineata, correct null model", xlab = "CORRECTED p-values")
dev.off()


diff_genes_1_1 <- res[which(res$padj < 0.05 & res$log2FoldChange  > 1),]
write.table(diff_genes_1_1, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_mod.txt", sep="\t")

diff_genes_1_2 <- res[which(res$padj < 0.05 & res$log2FoldChange < 0),]
write.table(diff_genes_1_2, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_mod.txt", sep="\t")


diff_genes_2_1 <- res_2[which(res_2$padj < 0.05 & res_2$log2FoldChange  > 1),]
write.table(diff_genes_2_1, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_mod_2.txt", sep="\t")

diff_genes_2_2 <- res_2[which(res_2$padj < 0.05 & res_2$log2FoldChange < 0),]
write.table(diff_genes_2_2, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_mod_2.txt", sep="\t")


diff_genes_G_H <- res_2G[which(res_2G$padj < 0.05 & res_2G$log2FoldChange  > 1),]
write.table(diff_genes_G_H, "/scr/bloodymary/rohit/Lacerta_viridis/DEseq/reps_high_DE_gonads.txt", sep="\t")

diff_genes_G_L <- res_2G[which(res_2G$padj < 0.05 & res_2G$log2FoldChange < 0),]
write.table(diff_genes_G_L, "/scr/bloodymary/rohit/Lacerta_viridis/DEseq/reps_low_DE_gonads.txt", sep="\t")

diff_genes_G_N <- res_2G[which(res_2G$padj < 0.05 & res_2G$log2FoldChange >= 0 & res_2G$log2FoldChange <= 1),]
write.table(diff_genes_G_L, "/scr/bloodymary/rohit/Lacerta_viridis/DEseq/reps_neutral_DE_gonads.txt", sep="\t")



diff_genes_g_G_H <- res_Group_Gonads[which(res_Group_Gonads$padj < 0.05 & res_Group_Gonads$log2FoldChange  > 1),]
write.table(diff_genes_G_H, "/scr/bloodymary/rohit/Lacerta_viridis/DEseq/reps_high_DE_group_gonads.txt", sep="\t")

diff_genes_g_G_L <- res_Group_Gonads[which(res_Group_Gonads$padj < 0.05 & res_Group_Gonads$log2FoldChange < 0),]
write.table(diff_genes_G_L, "/scr/bloodymary/rohit/Lacerta_viridis/DEseq/reps_low_DE_group_gonads.txt", sep="\t")

diff_genes_g_G_N <- res_Group_Gonads[which(res_Group_Gonads$padj < 0.05 & res_Group_Gonads$log2FoldChange >= 0 & res_Group_Gonads$log2FoldChange <= 1),]
write.table(diff_genes_G_L, "/scr/bloodymary/rohit/Lacerta_viridis/DEseq/reps_neutral_DE_group_gonads.txt", sep="\t")


ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)
summary(resLRT)
# adjusted p-value < 0.1
# LFC > 0 (up)     : 16079, 49% 
# LFC < 0 (down)   : 15722, 48% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# adjusted p-value < 0.1
# LFC > 0 (up)     : 11328, 34% 
# LFC < 0 (down)   : 11439, 34% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
resLRT<-resLRT[order(res$padj),]
head(resLRT)

DESeq2::plotMA(ddsLRT,ylim=c(-11,11),main='DESeq2 with Groups and Conditions - LRT')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_mod_LRT.png')
dev.off()

diff_genes_LRT <- resLRT[which(resLRT$padj < 0.05 & resLRT$log2FoldChange > 1),]
write.table(diff_genes_LRT, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_mod_LRT.txt", sep="\t")

diff_genes_2_LRT <- resLRT[which(resLRT$padj < 0.05 & resLRT$log2FoldChange < 0),]
write.table(diff_genes_2_LRT, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_mod_LRT.txt", sep="\t")


# res.B.only <- results(dds, contrast=c("group","Viridis","Bilineata"))
# res.B.only <-res.B.only[order(res.B.only$padj),]
# 
# res.H.only <- results(dds, contrast=c("condition","Viridis","Bilineata"))
# res.H.only <-res.B.only[order(res.B.only$padj),]




####################################################################################################################
### For Groups only
directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts",pattern="^new_features_")

sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),15))
sampleCondition = c(rep.int(c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2)),3))
sampleReplicate <- c(rep("1",10), rep("2",10),rep("3",10))

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,
                        group=sampleGroup, replicate=sampleReplicate)

#sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, group=sampleGroup)

ddsHTSeq_group<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group)
colData(ddsHTSeq_group)$group<-factor(colData(ddsHTSeq_group)$group, levels=c('Bilineata','Viridis'))

ddsHTSeq_group <- ddsHTSeq_group[, ddsHTSeq_group$group %in% c('Bilineata','Viridis') ]
ddsHTSeq_group$group <- relevel(ddsHTSeq_group$group, "Viridis")

dds_group <- estimateSizeFactors(ddsHTSeq_group)
idx_group <- rowSums( counts(dds_group, normalized=TRUE) >= 5 ) >= 3
dds_group <- dds_group[idx_group,]
dds_group <- DESeq(dds_group, minReplicatesForReplace=3)
#dds_group <- DESeq(ddsHTSeq_group)

## To check if normalisation worked with multidensity and multiecdf plots
GeneCounts_group <- counts(ddsHTSeq_group)
idx.nz_group <- apply(GeneCounts_group, 1, function(x) { all(x > 0)})
sum(idx.nz_group)
#15425

multidensity( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,],
              xlab="mean counts", xlim=c(0, 1000))
#multidensity( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,], 
#              xlab="mean counts", xlim=c(0, 1000), legend=("topright"))
dev.copy(png,'~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_multidensity.png')
dev.off()
multiecdf( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,],
              xlab="mean counts", xlim=c(0, 1000))
#multiecdf( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,], 
#              xlab="mean counts", xlim=c(0, 1000), legend=("topright"))
dev.copy(png,'~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_multiecdf.png')
dev.off()

## Pairwise comparison of samples with MD plot for mean and difference
pdf("~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_MDplot_pairwise.pdf")
MA.idx_group = t(combn(1:30, 2))
for( i in  seq_along( MA.idx_group[,1])){ 
  MDPlot(counts(ddsHTSeq_group, normalized = T)[idx.nz_group ,], 
         c(MA.idx_group[i,1],MA.idx_group[i,2]), 
         main = paste( colnames(ddsHTSeq_group)[MA.idx_group[i,1]], " vs ",
                       colnames(ddsHTSeq_group)[MA.idx_group[i,2]] ), ylim = c(-3,3))
}
dev.off()


countdata_group<-assay(ddsHTSeq_group)
coldata_group <- colData(ddsHTSeq_group)

nrow(ddsHTSeq_group)
#33824
ddsHTSeq_group <- ddsHTSeq_group[ rowSums(counts(ddsHTSeq_group)) > 1, ]
#33780
rld_group <- rlog(ddsHTSeq_group, blind=TRUE)
head(assay(rld_group), 3)

png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_scatterplot.png",width=1200,height=1100)
par( mfrow = c( 1, 2 ) )
ddsHTSeq_group <- estimateSizeFactors(ddsHTSeq_group)
plot(log2(counts(ddsHTSeq_group, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.3)
plot(assay(rld_group)[,1:2], pch=16, cex=0.3)
dev.off()

sampleDists_group <- dist( t( assay(rld_group) ) )
sampleDists_group
sampleDistMatrix_group <- as.matrix( sampleDists_group )
rownames(sampleDistMatrix_group) <- sub("*\\.txt", "", rownames(sampleDistMatrix_group))
rownames(sampleDistMatrix_group) <- sub(".*_rep", "rep", rownames(sampleDistMatrix_group))
colnames(sampleDistMatrix_group) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_heatplot_distmat.png",width=1200,height=1100)
pheatmap(sampleDistMatrix_group, clustering_distance_rows=sampleDists_group, clustering_distance_cols=sampleDists_group, col=colors)
dev.off()

poisd_group <- PoissonDistance(t(counts(ddsHTSeq_group)))
samplePoisDistMatrix_group <- as.matrix( poisd_group$dd )
rownames(samplePoisDistMatrix_group) <- sub("*\\.txt", "", rownames(samplePoisDistMatrix_group))
rownames(samplePoisDistMatrix_group) <- sub(".*_rep", "rep", rownames(samplePoisDistMatrix_group))
rownames(samplePoisDistMatrix_group) <- rownames(sampleDistMatrix_group)
colnames(samplePoisDistMatrix_group) <- NULL
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_heatplot_poisdist.png",width=1200,height=1100)
pheatmap(samplePoisDistMatrix_group, clustering_distance_rows=poisd_group$dd, clustering_distance_cols=poisd_group$dd, col=colors)
dev.off()


# ntop = 500
# Pvars_group <- rowVars(assay(rld_group))
# select_group <- order(Pvars_group, decreasing = TRUE)[seq_len(min(ntop, length(Pvars_group)))]
# PCA_group <- prcomp(t(assay(rld_group)[select, ]), scale = F)
# percentVar_group <- round(100*PCA_group$sdev^2/sum(PCA_group$sdev^2),1)


## To detect outliers
plotPCA(rld_group, intgroup = c("group","condition"))

data <- plotPCA(rld_group, intgroup = c("group","condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_PCA.png",width=1200,height=1100)
ggplot(data, aes(PC1, PC2, color=group.1, shape=condition)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()
dev.off()


outliers_group <- as.character(subset(colnames(ddsHTSeq_group), data$PC1 > 0))
outliers_group
#ddsHTSeq_group <- ddsHTSeq_group[, !(colnames(ddsHTSeq_group) %in% outliers)] 


# ddsHTSeq_group<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group*condition)
# rld_group <- rlog(ddsHTSeq_group, blind=TRUE)
# ggplot(data, aes(PC1, PC2, color=group.1, shape=condition)) + geom_point(size=3) + 
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()


png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_dispest.png",width=1200,height=1100)
plotDispEsts(estimateDispersions(ddsHTSeq_group))
dev.off()

#resultsNames(dds_group)
res_group<-results(dds_group)
res_group<-res_group[order(res_group$padj),]
head(res_group)
summary(res_group)
# out of 32603 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 2739, 8.4% 
# LFC < 0 (down)   : 3693, 11% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 7, 0.021% 


## Check for filtered-ones, pvalue and correct them accordingly
ddsHTSeq_group_2 <-  nbinomWaldTest(estimateDispersions(ddsHTSeq_group))
res_group_2 <- results(ddsHTSeq_group_2, pAdjustMethod = "BH")
table(res_group_2$padj < 0.1)
# FALSE  TRUE 
# 25842  5607 

png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_filterrej.png",width=1200,height=1100)
par(mfrow=c(1,2))
plot(metadata(res_group_2)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")
hist(res_group_2$pvalue, col = "lavender", main = "Viridis vs Bilineata", xlab = "p-values")
dev.off()


res_group_2 <- res_group_2[ !is.na(res_group_2$padj), ]
res_group_2 <- res_group_2[ !is.na(res_group_2$pvalue), ]
res_group_2 <- res_group_2[, -which(names(res_group_2) == "padj")]
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_pvaluecheck.png",width=1200,height=1100)
FDR.res_group_2 <- fdrtool(res_group_2$stat, statistic= "normal", plot = T)
dev.off()
FDR.res_group_2$param[1, "sd"]
#1.452559 

res_group_2[,"padj"]  <- p.adjust(FDR.res_group_2$pval, method = "BH")
hist(FDR.res_group_2$pval, col = "royalblue4", 
     main = "Viridis vs Bilineata, correct null model", xlab = "CORRECTED p-values")

table(res_group_2$padj < 0.1)
# FALSE  TRUE 
# 29960  1489 
plotMA(res_group_2, ylim=c(-6,6))

plotMA(dds_group,ylim=c(-5,5),main='DESeq2 with Groups only')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_group.png')
dev.off()

diff_genes_group <- res_group[which(res_group$padj < 0.05 & res_group$log2FoldChange > 1),]
write.table(diff_genes_group, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_group.txt", sep="\t")

diff_genes_group_2 <- res_group[which(res_group$padj < 0.05 & res_group$log2FoldChange < 0),]
write.table(diff_genes_group_2, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_group.txt", sep="\t")



######################################################################################################################
### For Conditions only
directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts",pattern="^new_features_")

sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),15))
sampleCondition = c(rep.int(c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2)),3))
sampleReplicate <- c(rep("1",10), rep("2",10),rep("3",10))

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,
                        group=sampleGroup, replicate=sampleReplicate)

#sampleCondition = c(rep.int(c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2)),3))
#sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

ddsHTSeq_cond<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq_cond)$condition<-factor(colData(ddsHTSeq_cond)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
ddsHTSeq_cond <- ddsHTSeq_cond[, ddsHTSeq_cond$condition %in% c('Brain','Gonads','Heart','Kidney','Liver') ]

ddsHTSeq_cond$condition <- relevel(ddsHTSeq_cond$condition, "Brain")

dds_cond <- estimateSizeFactors(ddsHTSeq_cond)
idx_cond <- rowSums( counts(dds_cond, normalized=TRUE) >= 5 ) >= 3
dds_cond <- dds_cond[idx_cond,]
dds_cond <- DESeq(dds_cond)

dds_cond<-DESeq(ddsHTSeq_cond)
res_cond<-results(dds_cond)
res_cond<-res_cond[order(res_cond$padj),]
head(res_cond)

dds_cond_LRT <- DESeq(ddsHTSeq_cond, test="LRT", reduced= ~ 1)
res_cond_LRT <- results(dds_cond_LRT)

DESeq2::plotMA(dds_cond,ylim=c(-11,11),main='DESeq2 with Conditions only')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_cond.png')
dev.off()

diff_genes_cond <- res_cond[which(res_cond$padj < 0.05 & res_cond$log2FoldChange > 1),]
write.table(diff_genes_cond, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_cond.txt", sep="\t")

diff_genes_cond_2 <- res_cond[which(res_cond$padj < 0.05 & res_cond$log2FoldChange < 0),]
write.table(diff_genes_cond_2, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_cond.txt", sep="\t")




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


png(file="~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_cond_pairwise.png",width=1100,height=800)
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
#dev.copy(png,'~/Downloads/scripts/reps_deseq2_MAplot_cond_pairwise.png')
dev.off()



#####################################################################################################################
# Groups with tissue interaction

directory <- file.path("/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts")
sampleFiles <- list.files(path="/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts",pattern="^new_features_")


### For Condition and Group
sampleGroup <- factor(rep.int(c("Bilineata","Viridis"),15))
sampleCondition = c(rep.int(c(rep("Brain",2), rep("Gonads",2),rep("Heart",2), rep("Kidney",2),rep("Liver",2)),3))
sampleReplicate <- c(rep("1",10), rep("2",10),rep("3",10))
#sampleReplicate <- factor(rep.int(c("1","2","3"),5))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,
                        group=sampleGroup, replicate=sampleReplicate)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design =~ group+condition)
#ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design =~ group+condition+group:condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
colData(ddsHTSeq)$group<-factor(colData(ddsHTSeq)$group, levels=c('Bilineata','Viridis'))
#colData(ddsHTSeq)$replicate<-factor(colData(ddsHTSeq)$replicate, levels=c('1','2','3'))

ddsHTSeq <- ddsHTSeq[, ddsHTSeq$group %in% c('Bilineata','Viridis') ]
ddsHTSeq <- ddsHTSeq[, ddsHTSeq$condition %in% c('Brain','Gonads','Heart','Kidney','Liver') ]
#ddsHTSeq <- ddsHTSeq[, ddsHTSeq$replicate %in% c('1','2','3') ]
ddsHTSeq$group <- relevel(ddsHTSeq$group, "Viridis")

nrow(ddsHTSeq)
#33824
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]
#33780
rld <- rlog(ddsHTSeq, blind=TRUE)
head(assay(rld), 3)


## To detect outliers
plotPCA(ddsHTSeq, intgroup = c("group","condition"))
data <- plotPCA(rld, intgroup = c("group","condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
png(file="~/Downloads/scripts/differential_expression/reps_deseq2_PCA_int.png",width=1200,height=1100)
ggplot(data, aes(PC1, PC2, color=group.1, shape=condition)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()
dev.off()

outliers <- as.character(subset(colnames(ddsHTSeq), data$PC1 > 0))
outliers

dds_int <- estimateSizeFactors(ddsHTSeq)
idx <- rowSums( counts(dds_int, normalized=TRUE) >= 10 ) >= 3
#idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds_int <- dds_int[idx,]
dds_int <- DESeq(dds_int)

GeneCounts <- counts(ddsHTSeq)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)
#15425


## To check if normalisation worked with multidensity and multiecdf plots, pairwise MDplots
multidensity( counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_multidensity_mod_int.png')
dev.off()
multiecdf( counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1000))
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_multiecdf_mod_int.png')
dev.off()

#png(file="~/Downloads/scripts/differential_expression/reps_deseq2_MDplot_pairwise_int.png",width=1200,height=1100)
pdf("~/Downloads/scripts/differential_expression/reps_deseq2_MDplot_pairwise_int.pdf")
MA.idx = t(combn(1:30, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,], 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(ddsHTSeq)[MA.idx[i,1]], " vs ",
                       colnames(ddsHTSeq)[MA.idx[i,2]] ), ylim = c(-3,3))
}
dev.off()


png(file="~/Downloads/scripts/differential_expression/reps_deseq2_dispest_int.png",width=1200,height=1100)
plotDispEsts(estimateDispersions(estimateSizeFactors(ddsHTSeq)))
dev.off()


resultsNames(dds_int)
#res_int<-results(dds_int, contrast=list(c("groupViridis","groupBilineata")))
res_int<-results(dds_int, name="groupViridis")
#res<-results(dds, list( c("Intercept") ))
summary(res_int)
### Interaction * with zero-count removal (>=10)
# out of 30053 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 7234, 24% 
# LFC < 0 (down)   : 8118, 27% 
# outliers [1]     : 46, 0.15% 
# low counts [2]   : 0, 0% 
# (mean count < 1)


res_int<-res_int[order(res_int$padj),]
head(res_int)


diff_genes_1_1_int <- res_int[which(res_int$padj < 0.05 & res_int$log2FoldChange  > 1),]
write.table(diff_genes_1_1_int, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_mod_int.txt", sep="\t")

diff_genes_1_2 <- res_int[which(res_int$padj < 0.05 & res_int$log2FoldChange < 0),]
write.table(diff_genes_1_2, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_mod_int.txt", sep="\t")


plotMA(dds_int,ylim=c(-11,11),main='DESeq2 with Groups and Conditions')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_mod_int.png')
dev.off()

resLFC1_int <- results(dds_int, lfcThreshold=1)
topGene <- rownames(resLFC1_int)[which.min(resLFC1_int$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="darkblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="darkblue")
})


## Check for filtered-ones, pvalue and correct them accordingly
ddsHTSeq_2 <-  nbinomWaldTest(estimateDispersions(estimateSizeFactors(ddsHTSeq)))
res_2_int <- results(ddsHTSeq_2, pAdjustMethod = "BH", name="groupViridis")
table(res_2$padj < 0.1)
# FALSE  TRUE 
# 30161 2309


png(file="~/Downloads/scripts/differential_expression/reps_deseq2_filterrej_int.png",width=1200,height=1100)
par(mfrow=c(1,2))
plot(metadata(res_2)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")
hist(res_2$pvalue, col = "lavender", main = "WT vs Deletion", xlab = "p-values")
dev.off()

res_2_int <- res_2_int[ !is.na(res_2_int$padj), ]
res_2_int <- res_2_int[ !is.na(res_2_int$pvalue), ]
res_2_int <- res_2_int[, -which(names(res_2_int) == "padj")]
png(file="~/Downloads/scripts/differential_expression/reps_deseq2_pvaluecheck_int.png",width=1200,height=1100)
FDR.res_2_int <- fdrtool(res_2_int$stat, statistic= "normal", plot = T)
dev.off()
FDR.res_2_int$param[1, "sd"]
# sd 
# 2.603573

png(file="~/Downloads/scripts/differential_expression/reps_deseq2_pvaluecorrected_int.png",width=1200,height=1100)
res_2_int[,"padj"]  <- p.adjust(FDR.res_2_int$pval, method = "BH")
hist(FDR.res_2_int$pval, col = "royalblue4", 
     main = "Viridis vs Bilineata, correct null model", xlab = "CORRECTED p-values")
dev.off()


diff_genes_2_1_int <- res_2_int[which(res_2_int$padj < 0.05 & res_2_int$log2FoldChange  > 1),]
write.table(diff_genes, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_mod_2_int.txt", sep="\t")

diff_genes_2_2_int <- res_2_int[which(res_2_int$padj < 0.05 & res_2_int$log2FoldChange < 0),]
write.table(diff_genes_2_2_int, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_mod_2_int.txt", sep="\t")


ddsLRT_int <- DESeq(dds_int, test="LRT", reduced= ~ 1)
resLRT_int <- results(ddsLRT_int, name="groupViridis")
summary(resLRT_int)
# out of 32610 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 14284, 44% 
# LFC < 0 (down)   : 17082, 52% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 1)

resLRT_int<-resLRT_int[order(res_int$padj),]
head(resLRT_int)

DESeq2::plotMA(ddsLRT_int,ylim=c(-11,11),main='DESeq2 with Groups and Conditions - LRT')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_mod_LRT_int.png')
dev.off()

diff_genes_LRT_int <- resLRT[which(resLRT$padj < 0.05 & resLRT_int$log2FoldChange > 1),]
write.table(diff_genes_LRT_int, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_differentially_expressed_mod_LRT_int.txt", sep="\t")

diff_genes_2_LRT_int <- resLRT[which(resLRT_int$padj < 0.05 & resLRT_int$log2FoldChange < 0),]
write.table(diff_genes_2_LRT_int, "/scr/k70san/rohit/Lacertidae/splicing/unpooled/counts/reps_lowly_expressed_mod_LRT_int.txt", sep="\t")


