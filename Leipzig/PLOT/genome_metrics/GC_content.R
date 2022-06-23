library(ggplot2)
library(gridExtra)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/GC-content/")

data1 <- read.table(file = "distribution_GC.txt",header = F,sep = "\t")

ggplot(data1, aes(x=data1$V1, fill=data1$V2)) + geom_density(alpha=.3) + scale_x_continuous(limits = c(0.35, 0.65)) +
  theme_bw() + theme(legend.title=element_blank()) #+ facet_wrap(~data1$V2)

hist(data1$V1, freq = T, xlab = 'Fraction of GC-content', xlim = c(0.3, 0.65), 
     ylab = 'Density distribution', main = 'Histogram of GC-content in L.viridis with Kernel Density Plot')
lines(density(data1$V1, na.rm = T, from = 0.35, to = 0.65))


GC_all <- subset(data1, data1$V2 == "ALL")
GC_notw <- subset(data1, data1$V2 == "Not-W")
GC_notz <- subset(data1, data1$V2 == "Not-Z")
GC_notzw <- subset(data1, data1$V2 == "Not-ZW")
GC_w <- subset(data1, data1$V2 == "W-only")
GC_z <- subset(data1, data1$V2 == "Z-only")
GC_zw <- subset(data1, data1$V2 == "ZW-both")

# GC_test <- ggplot(GC_all, aes(x=GC_all$V1, fill=GC_all$V2)) + geom_density(alpha=.2) + 
#   scale_x_continuous(limits = c(0.35, 0.65)) + theme_bw() + theme(legend.title=element_blank())
# 
# GC_1 <- ggplot(GC_all, aes(x=GC_all$V1, fill=GC_all$V2)) + geom_density(alpha=.1) + 
#   geom_histogram(binwidth =.0015, alpha=.3) +
#   scale_x_continuous(limits = c(0.35, 0.65)) + theme_bw() + theme(legend.title=element_blank())


GC_test <- ggplot(GC_all, aes(x=GC_all$V1, fill=GC_all$V2)) + geom_density(alpha=.2) + 
  theme_bw() + theme(legend.title=element_blank())

GC_1 <- ggplot(GC_all, aes(x=GC_all$V1, fill=GC_all$V2)) + geom_density(alpha=.1) + 
  geom_histogram(binwidth =.0015, alpha=.3) +
  theme_bw() + theme(legend.title=element_blank())

GC_2 <- ggplot(GC_notw, aes(x=GC_notw$V1, fill=GC_notw$V2)) + geom_density(alpha=.1) + 
  geom_histogram(binwidth =.0015, alpha=.4) +
  theme_bw() + theme(legend.title=element_blank())

GC_3 <- ggplot(GC_notz, aes(x=GC_notz$V1, fill=GC_notz$V2)) + geom_density(alpha=.1) + 
  geom_histogram(binwidth =.0015, alpha=.4) +
  theme_bw() + theme(legend.title=element_blank())

GC_4 <- ggplot(GC_notzw, aes(x=GC_notzw$V1, fill=GC_notzw$V2)) + geom_density(alpha=.1) + 
  geom_histogram(binwidth =.0015, alpha=.5) +
  theme_bw() + theme(legend.title=element_blank())

GC_5 <- ggplot(GC_w, aes(x=GC_w$V1, fill=GC_w$V2)) + geom_density(alpha=.1) + geom_histogram(binwidth =.0015, alpha=.6) +
  theme_bw() + theme(legend.title=element_blank())

GC_6 <- ggplot(GC_z, aes(x=GC_z$V1, fill=GC_z$V2)) + geom_density(alpha=.1) + geom_histogram(binwidth =.0015, alpha=.7) +
  theme_bw() + theme(legend.title=element_blank())

GC_7 <- ggplot(GC_zw, aes(x=GC_zw$V1, fill=GC_zw$V2)) + geom_density(alpha=.1) + geom_histogram(binwidth =.0015, alpha=.8) +
  theme_bw() + theme(legend.title=element_blank())


grid.arrange(GC_test,GC_1, GC_2, GC_3, GC_4, GC_5, GC_6, GC_7, ncol=2)
