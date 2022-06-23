
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

setwd("/scr/k61san/nowicklab/Lacerta/DEseq")
# Load directory and get mapping-files


#### Working version

directory <- file.path("/scr/k61san/nowicklab/Lacerta/DEseq")
sampleFiles <- list.files(path="/scr/k61san/nowicklab/Lacerta/DEseq",pattern="^new_features_")


### For Condition and Group
sampleGroup <- factor(rep.int(c(rep("Adriatic",4), rep("Viridis",3)),5))
sampleCondition <- c(c(rep("Brain",7), rep("Gonads",7),rep("Heart",7), rep("Kidney",7),rep("Liver",7)))
sampleReplicate <- factor(rep.int(c(rep(c("1","2"),3), c("3")), 5))
sampleSex <- factor(rep.int(c(rep("Female",2), rep("Male",2), rep("Female",3)), 5))

#sampleReplicate <- factor(rep.int(c("1","2","3"),5))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,
                        group=sampleGroup, replicate=sampleReplicate, sex=sampleSex)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design =~ group*condition)
#ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design =~ group+condition+group:condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('Brain','Gonads','Heart','Kidney','Liver'))
colData(ddsHTSeq)$group<-factor(colData(ddsHTSeq)$group, levels=c('Adriatic','Viridis'))
#colData(ddsHTSeq)$replicate<-factor(colData(ddsHTSeq)$replicate, levels=c('1','2','3'))

ddsHTSeq <- ddsHTSeq[, ddsHTSeq$group %in% c('Adriatic','Viridis') ]
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
png(file="~/Downloads/scripts/differential_expression/reps_deseq2_PCA_virAdr.png",width=1200,height=1100)
ggplot(data, aes(PC1, PC2, color=group.1, shape=condition)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()
dev.off()

outliers <- as.character(subset(colnames(ddsHTSeq), data$PC1 > 0))
outliers

dds <- estimateSizeFactors(ddsHTSeq)
idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 3
#idx <- rowSums( counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

GeneCounts <- counts(ddsHTSeq)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)
#13162


## To check if normalisation worked with multidensity and multiecdf plots, pairwise MDplots
multidensity( counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))
dev.copy(png,'~/Downloads/scripts/reps_deseq2_multidensity_mod_virAdr.png')
dev.off()
multiecdf( counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1000))
dev.copy(png,'~/Downloads/scripts/reps_deseq2_multiecdf_mod_virAdr.png')
dev.off()

#png(file="~/Downloads/scripts/differential_expression/reps_deseq2_MDplot_pairwise.png",width=1200,height=1100)
pdf("~/Downloads/scripts/differential_expression/reps_deseq2_MDplot_pairwise_virAdr.pdf")
MA.idx = t(combn(1:30, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(estimateSizeFactors(ddsHTSeq), normalized = T)[idx.nz ,], 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(ddsHTSeq)[MA.idx[i,1]], " vs ",
                       colnames(ddsHTSeq)[MA.idx[i,2]] ), ylim = c(-3,3))
}
dev.off()


png(file="~/Downloads/scripts/differential_expression/reps_deseq2_dispest_virAdr.png",width=1200,height=1100)
plotDispEsts(estimateDispersions(estimateSizeFactors(ddsHTSeq)))
dev.off()


# idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
# dds <- dds[idx,]
# dds <- DESeq(dds)

res<-results(dds)
summary(res)
### Interaction * with zero-count removal (>=5)
# out of 33319 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 1467, 4.4% 
# LFC < 0 (down)   : 1256, 3.8% 
# outliers [1]     : 94, 0.28% 
# low counts [2]   : 2571, 7.7% 
# (mean count < 7)
### Interaction * with zero-count removal (>=10)
# out of 32599 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 1465, 4.5% 
# LFC < 0 (down)   : 1256, 3.9% 
# outliers [1]     : 87, 0.27% 
# low counts [2]   : 1888, 5.8% 
# (mean count < 7)
### Interaction * with zero-count removal (>=100)
# out of 20590 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 1157, 5.6% 
# LFC < 0 (down)   : 1089, 5.3% 
# outliers [1]     : 42, 0.2% 
# low counts [2]   : 0, 0% 
# (mean count < 13)


res<-res[order(res$padj),]
head(res)


plotMA(dds,ylim=c(-11,11),main='DESeq2 with Groups and Conditions')
dev.copy(png,'~/Downloads/scripts/reps_deseq2_MAplot_mod_virAdr.png')
dev.off()

resLFC1 <- results(dds, lfcThreshold=1)
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="darkblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="darkblue")
})


## Check for filtered-ones, pvalue and correct them accordingly
ddsHTSeq_2 <-  nbinomWaldTest(estimateDispersions(estimateSizeFactors(ddsHTSeq)))
res_2 <- results(ddsHTSeq_2, pAdjustMethod = "BH")
table(res_2$padj < 0.05)
# FALSE  TRUE 
# 28698  1751


png(file="~/Downloads/scripts/differential_expression/reps_deseq2_filterrej_virAdr.png",width=1200,height=1100)
par(mfrow=c(1,2))
plot(metadata(res_2)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")
hist(res_2$pvalue, col = "lavender", main = "WT vs Deletion", xlab = "p-values")
dev.off()

res_2 <- res_2[ !is.na(res_2$padj), ]
res_2 <- res_2[ !is.na(res_2$pvalue), ]
res_2 <- res_2[, -which(names(res_2) == "padj")]
png(file="~/Downloads/scripts/differential_expression/reps_deseq2_pvaluecheck_virAdr.png",width=1200,height=1100)
FDR.res_2 <- fdrtool(res_2$stat, statistic= "normal", plot = T)
dev.off()
FDR.res_2$param[1, "sd"]
# sd 
# 1.443763

png(file="~/Downloads/scripts/differential_expression/reps_deseq2_pvaluecorrected_virAdr.png",width=1200,height=1100)
res_2[,"padj"]  <- p.adjust(FDR.res_2$pval, method = "BH")
hist(FDR.res_2$pval, col = "royalblue4", 
     main = "Viridis vs Adriatic, correct null model", xlab = "CORRECTED p-values")
dev.off()


diff_genes <- res[which(res$padj < 0.05 & res$log2FoldChange  > 1),]
write.table(diff_genes, "reps_differentially_expressed_mod_virAdr.txt", sep="\t")

diff_genes_2 <- res[which(res$padj < 0.05 & res$log2FoldChange < 0),]
write.table(diff_genes_2, "reps_lowly_expressed_mod_virAdr.txt", sep="\t")


diff_genes_2_1 <- res_2[which(res_2$padj < 0.05 & res_2$log2FoldChange  > 1),]
write.table(diff_genes, "reps_differentially_expressed_mod_2_virAdr.txt", sep="\t")

diff_genes_2_2 <- res_2[which(res_2$padj < 0.05 & res_2$log2FoldChange < 0),]
write.table(diff_genes_2, "reps_lowly_expressed_mod_2_virAdr.txt", sep="\t")


ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)
summary(resLRT)
# out of 20590 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 9218, 45% 
# LFC < 0 (down)   : 10762, 52% 
# outliers [1]     : 42, 0.2% 
# low counts [2]   : 0, 0% 
# (mean count < 13)


resLRT<-resLRT[order(res$padj),]
head(resLRT)

DESeq2::plotMA(ddsLRT,ylim=c(-11,11),main='DESeq2 with Groups and Conditions - LRT')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_mod_LRT_virAdr.png')
dev.off()

diff_genes_LRT <- resLRT[which(resLRT$padj < 0.05 & resLRT$log2FoldChange > 1),]
write.table(diff_genes_LRT, "reps_differentially_expressed_mod_LRT_virAdr.txt", sep="\t")

diff_genes_2_LRT <- resLRT[which(resLRT$padj < 0.05 & resLRT$log2FoldChange < 0),]
write.table(diff_genes_2_LRT, "reps_lowly_expressed_mod_LRT_virAdr.txt", sep="\t")


####################################################################################################################
### For Groups only
directory <- file.path("/scr/k61san/nowicklab/Lacerta/DEseq")
sampleFiles <- list.files(path="/scr/k61san/nowicklab/Lacerta/DEseq",pattern="^new_features_")

sampleGroup <- factor(rep.int(c(rep("Adriatic",4), rep("Viridis",3)),5))
sampleCondition <- c(c(rep("Brain",7), rep("Gonads",7),rep("Heart",7), rep("Kidney",7),rep("Liver",7)))
sampleReplicate <- factor(rep.int(c(rep(c("1","2"),3), c("3")), 5))
sampleSex <- factor(rep.int(c(rep("Female",2), rep("Male",2), rep("Female",3)), 5))

#sampleReplicate <- factor(rep.int(c("1","2","3"),5))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,
                        group=sampleGroup, replicate=sampleReplicate, sex=sampleSex)

#sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, group=sampleGroup)

ddsHTSeq_group<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~group)
colData(ddsHTSeq_group)$group<-factor(colData(ddsHTSeq_group)$group, levels=c('Adriatic','Viridis'))

ddsHTSeq_group <- ddsHTSeq_group[, ddsHTSeq_group$group %in% c('Adriatic','Viridis') ]
ddsHTSeq_group$group <- relevel(ddsHTSeq_group$group, "Viridis")

dds_group <- estimateSizeFactors(ddsHTSeq_group)
idx_group <- rowSums( counts(dds_group, normalized=TRUE) >= 5 ) >= 3
dds_group <- dds_group[idx_group,]
dds_group <- DESeq(dds_group, minReplicatesForReplace=3)


## To check if normalisation worked with multidensity and multiecdf plots
GeneCounts_group <- counts(ddsHTSeq_group)
idx.nz_group <- apply(GeneCounts_group, 1, function(x) { all(x > 0)})
sum(idx.nz_group)
#13162

multidensity( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,],
              xlab="mean counts", xlim=c(0, 1000))
#multidensity( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,], 
#              xlab="mean counts", xlim=c(0, 1000), legend=("topright"))
dev.copy(png,'~/Downloads/scripts/new_deseq2_replicated_groups_multidensity_virAdr.png')
dev.off()
multiecdf( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,],
              xlab="mean counts", xlim=c(0, 1000))
#multiecdf( counts(estimateSizeFactors(ddsHTSeq_group), normalized = T)[idx.nz_group ,], 
#              xlab="mean counts", xlim=c(0, 1000), legend=("topright"))
dev.copy(png,'~/Downloads/scripts/new_deseq2_replicated_groups_multiecdf_virAdr.png')
dev.off()

## Pairwise comparison of samples with MD plot for mean and difference
# pdf("~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_MDplot_pairwise_virAdr.pdf")
# MA.idx_group = t(combn(1:30, 2))
# for( i in  seq_along( MA.idx_group[,1])){ 
#   MDPlot(counts(ddsHTSeq_group, normalized = T)[idx.nz_group ,], 
#          c(MA.idx_group[i,1],MA.idx_group[i,2]), 
#          main = paste( colnames(ddsHTSeq_group)[MA.idx_group[i,1]], " vs ",
#                        colnames(ddsHTSeq_group)[MA.idx_group[i,2]] ), ylim = c(-3,3))
# }
# dev.off()


countdata_group<-assay(ddsHTSeq_group)
coldata_group <- colData(ddsHTSeq_group)

nrow(ddsHTSeq_group)
#33824
ddsHTSeq_group <- ddsHTSeq_group[ rowSums(counts(ddsHTSeq_group)) > 1, ]
#33809
rld_group <- rlog(ddsHTSeq_group, blind=TRUE)
head(assay(rld_group), 3)

png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_scatterplot_virAdr.png",width=1200,height=1100)
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
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_heatplot_distmat_virAdr.png",width=1200,height=1100)
pheatmap(sampleDistMatrix_group, clustering_distance_rows=sampleDists_group, clustering_distance_cols=sampleDists_group, col=colors)
dev.off()

poisd_group <- PoissonDistance(t(counts(ddsHTSeq_group)))
samplePoisDistMatrix_group <- as.matrix( poisd_group$dd )
rownames(samplePoisDistMatrix_group) <- sub("*\\.txt", "", rownames(samplePoisDistMatrix_group))
rownames(samplePoisDistMatrix_group) <- sub(".*_rep", "rep", rownames(samplePoisDistMatrix_group))
rownames(samplePoisDistMatrix_group) <- rownames(sampleDistMatrix_group)
colnames(samplePoisDistMatrix_group) <- NULL
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_heatplot_poisdist_virAdr.png",width=1200,height=1100)
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
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_PCA_virAdr.png",width=1200,height=1100)
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


png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_dispest_virAdr.png",width=1200,height=1100)
plotDispEsts(estimateDispersions(ddsHTSeq_group))
dev.off()


res_group<-results(dds_group)
res_group<-res_group[order(res_group$padj),]
head(res_group)
summary(res_group)
# out of 33319 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 2289, 6.9% 
# LFC < 0 (down)   : 4165, 13% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 


## Check for filtered-ones, pvalue and correct them accordingly
ddsHTSeq_group_2 <-  nbinomWaldTest(estimateDispersions(ddsHTSeq_group))
res_group_2 <- results(ddsHTSeq_group_2, pAdjustMethod = "BH")
table(res_group_2$padj < 0.1)
# FALSE  TRUE 
# 24154  5476 

png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_filterrej_virAdr.png",width=1200,height=1100)
par(mfrow=c(1,2))
plot(metadata(res_group_2)$filterNumRej, type="b", xlab="quantiles of 'baseMean'",
     ylab="number of rejections")
hist(res_group_2$pvalue, col = "lavender", main = "Viridis vs Adriatic", xlab = "p-values")
dev.off()


res_group_2 <- res_group_2[ !is.na(res_group_2$padj), ]
res_group_2 <- res_group_2[ !is.na(res_group_2$pvalue), ]
res_group_2 <- res_group_2[, -which(names(res_group_2) == "padj")]
png(file="~/Downloads/scripts/differential_expression/new_deseq2_replicated_groups_pvaluecheck_virAdr.png",width=1200,height=1100)
FDR.res_group_2 <- fdrtool(res_group_2$stat, statistic= "normal", plot = T)
dev.off()
FDR.res_group_2$param[1, "sd"]
#1.63288

res_group_2[,"padj"]  <- p.adjust(FDR.res_group_2$pval, method = "BH")
hist(FDR.res_group_2$pval, col = "royalblue4", 
     main = "Viridis vs Adriatic, correct null model", xlab = "CORRECTED p-values")

table(res_group_2$padj < 0.05)
# FALSE  TRUE 
# 29188  442
plotMA(res_group_2, ylim=c(-6,6))

plotMA(dds_group,ylim=c(-5,5),main='DESeq2 with Groups only')
dev.copy(png,'~/Downloads/scripts/differential_expression/reps_deseq2_MAplot_group_virAdr.png')
dev.off()

diff_genes_group <- res_group[which(res_group$padj < 0.05 & res_group$log2FoldChange > 1),]
write.table(diff_genes_group, "reps_differentially_expressed_group_virAdr.txt", sep="\t")

diff_genes_group_2 <- res_group[which(res_group$padj < 0.05 & res_group$log2FoldChange < 0),]
write.table(diff_genes_group_2, "reps_lowly_expressed_group_virAdr.txt", sep="\t")



######################################################################################################################
### For Conditions only
directory <- file.path("/scr/k61san/nowicklab/Lacerta/DEseq")
sampleFiles <- list.files(path="/scr/k61san/nowicklab/Lacerta/DEseq",pattern="^new_features_")

sampleGroup <- factor(rep.int(c(rep("Adriatic",4), rep("Viridis",3)),5))
sampleCondition <- c(c(rep("Brain",7), rep("Gonads",7),rep("Heart",7), rep("Kidney",7),rep("Liver",7)))
sampleReplicate <- factor(rep.int(c(rep(c("1","2"),3), c("3")), 5))
sampleSex <- factor(rep.int(c(rep("Female",2), rep("Male",2), rep("Female",3)), 5))

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition,
                        group=sampleGroup, replicate=sampleReplicate, sex=sampleSex)


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
dev.copy(png,'~/Downloads/scripts/reps_deseq2_MAplot_cond_virAdr.png')
dev.off()

diff_genes_cond <- res_cond[which(res_cond$padj < 0.05 & res_cond$log2FoldChange > 1),]
write.table(diff_genes_cond, "reps_differentially_expressed_cond_virAdr.txt", sep="\t")

diff_genes_cond_2 <- res_cond[which(res_cond$padj < 0.05 & res_cond$log2FoldChange < 0),]
write.table(diff_genes_cond_2, "reps_lowly_expressed_cond_virAdr.txt", sep="\t")



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


png(file="~/Downloads/scripts/reps_deseq2_MAplot_cond_pairwise_virAdr.png",width=1100,height=800)
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
dev.off()


