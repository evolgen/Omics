
library(ggplot2)
library(gridExtra)
library(reshape)
require(grid)
library(ggpubr)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/Selection")

dt_MA <- read.table(file = "divergence_time_MA.txt",header = F,sep = "\t")
dt_GMYN <- read.table(file = "divergence_time_GMYN.txt",header = F,sep = "\t")

#plot1 <- plot(density(dt_MA$V1), xlim = c(1e+5, 5e+6))
#plot2 <- plot(density(dt_GMYN$V1), xlim = c(1e+5, 5e+6))

plot1 <- ggplot(dt_MA, aes(dt_MA$V1)) + geom_density() + xlim(1e+4, 5e+6)
plot2 <- ggplot(dt_GMYN, aes(dt_GMYN$V1)) + geom_density() + xlim(1e+4, 5e+6)

grid.arrange(plot1, plot2, ncol=2)
