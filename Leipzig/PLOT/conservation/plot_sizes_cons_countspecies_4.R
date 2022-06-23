setwd("~/Downloads/scripts/conservation/")

library(ggplot2)
library(gridExtra)

cons_xen<-read.table("./sizes_xen_count_species_4", header=FALSE)
cons_gal<-read.table("./sizes_gal_count_species_4", header=FALSE)


summary(cons_xen$V1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   34.00   51.00  141.00   97.75 5049.00 
sd(cons_xen$V1)
# 324.4312


summary(cons_gal$V1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   22.75   92.00  132.10  203.20  830.00 
sd(cons_gal$V1)
# 136.1417


plot1 <- ggplot(cons_xen,aes(x=cons_xen$V1)) + geom_histogram(binwidth = 10)
plot2 <- ggplot(cons_gal,aes(x=cons_gal$V1)) + geom_histogram(binwidth = 10)

grid.arrange(plot1, plot2, ncol=1)

