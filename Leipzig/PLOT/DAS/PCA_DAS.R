
# Load data
data(iris)
head(iris, 3)

# log transform 
log.ir <- log(iris[, 1:4])
ir.species <- iris[, 5]

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
ir.pca <- prcomp(log.ir, center = TRUE, scale. = TRUE)

# print method
print(ir.pca)

# plot method
plot(ir.pca, type = "l")

# summary method
summary(ir.pca)

# Predict PCs
predict(ir.pca, newdata=tail(log.ir, 2))


library(devtools)
install_github("vqv/ggbiplot", force = TRUE)
library(ggbiplot)

g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, groups = ir.species, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)

library(data.table)
library(reshape2)
data1 <- fread("/scr/k70san/rohit/Lacertidae/splicing/pooled/ALL/HTseq/TRIAD/AlternativeSplicing/VIR_diego_species_three/junction_table.txt")
head(data1,3)
data1_1 <- data1[,-c(2,44)]
data1_1[data1_1 == 0] <- NA
data1_2 <- na.omit(data1_1)

data1_3 <- melt(data1_2,id=c("junction","geneID"))
log.data <- log(data1_3[, 4])
data.species <- unique(data1_2$variable)
data.pca <- prcomp(log.data, center = TRUE, scale. = TRUE)
print(data.pca)

plot(data.pca, type = "l")
summary(data.pca)
predict(data.pca, newdata=tail(log.data, 2))


log.data.2 <- log(data1_2[, 2:41])
data.pca <- prcomp(log.data.2, center = TRUE, scale. = TRUE)
print(data.pca)

plot(data.pca, type = "l")
summary(data.pca)
predict(data.pca, newdata=tail(log.data.2, 2))


