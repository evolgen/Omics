setwd("/scr/k61san/nowicklab/Lacerta/self_vir")

library(pastecs)
library(ggplot2)
#remove.packages("ggplot2")
library(gridExtra)
#library(grid)
library(RColorBrewer)
library(nortest)


data_SegDupl<-read.table("/scr/k61san/nowicklab/Lacerta/self_vir/SEGM_DUP.lacerta.txt", 
                        header=T, sep = "\t")
stat.desc(data_SegDupl$viridis)
stat.desc(data_SegDupl$bilineata)
var.test(data_SegDupl$viridis, data_SegDupl$bilineata)
# F = 1.1318, num df = 327, denom df = 191, p-value = 0.3453
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.8750475 1.4515531
# sample estimates:
#   ratio of variances 
# 1.13183 

wilcox.test(data_SegDupl$viridis, data_SegDupl$bilineata, paired=FALSE, alternative = c("two.sided"))
# W = 32474, p-value = 0.551
# alternative hypothesis: true location shift is not equal to 0

