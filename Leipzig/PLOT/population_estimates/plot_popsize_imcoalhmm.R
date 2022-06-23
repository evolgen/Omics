setwd("~/Downloads/scripts/population_estimates/")
library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)

# 'isolation.period', 'migration.period','theta', 'rho', 'migration', 'log.likelihood'

#ImchDat_ts<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/I_M_times_scaled", header=FALSE)

dat<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/all_mod_models_I_M", header=FALSE)

ggplot(dat) + geom_density(aes(x=dat$V1), fill="red", alpha=0.1) + geom_density(aes(x=dat$V2), fill = "lightblue", alpha=0.4) +
  xlim(-0.001,0.005) +theme_bw()


dat2<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/I_M_times_scaled", header=FALSE)

ggplot(dat2) + geom_density(aes(x=dat2$V1), fill="red", alpha=0.1, show.legend=TRUE) + 
  geom_density(aes(x=dat2$V2), fill = "lightblue", alpha=0.4) +
  theme_bw() + scale_x_continuous(breaks=c(4000,30000,75000,200000,450000), limits = c(-50000,450000) ) +
  ylab("Frequency") + xlab("Time in Years") + theme_bw()


dat3 <- gather(dat2)
bitmap("/homes/biertank/rohit/Downloads/scripts/population_estimates/IMcoalHMM_estimates.tiff", height = 8, width = 12, units = 'in', res=1600)

ggplot(dat3, aes(x=value, fill=factor(key,labels=c("Isolation-time","Migration-time")), 
  alpha=0.2)) + geom_density() +
  scale_x_continuous(breaks=c(4000,30000,75000,200000,450000), limits = c(-50000,450000)) + theme_bw() +
  guides(fill=guide_legend(title="Parameter-estimates"), alpha=FALSE) + 
  labs(title="Estimation of Split-times in Lacertids with Isolation-with-Migration model of IMCoalHMM") +
  ylab("Frequency") + xlab("Time in Years")# + guides(alpha = FALSE) 
                                  
dev.copy(png,"/homes/biertank/rohit/Downloads/scripts/population_estimates/IMcoalHMM_estimates.png", height = 8, width = 12, units = 'in', res=1600)
dev.off()



#m1 + geom_freqpoly(fill=c("light-blue")) + xlim(0,1.5e+12) + ylim(0,3500)
m1 + geom_area(aes(y = ..count..), fill="lightblue", stat = "bin") + xlim(-0.5,1.5e+12) + ylim(0,3500) +
  geom_smooth(method="gam")
m2 + geom_area(aes(y = ..count..), fill="red", stat = "bin") + xlim(-0.5,1.5e+12) + ylim(0,3500) +
  geom_smooth(method="gam")

# 'isolation.period', 'migration.period','theta', 'rho', 'migration', 'log.likelihood'

# setwd("/scr/k61san/nowicklab/Lacerta/IMcoal/")
# library(ggplot2)
# ImchDat_I_only<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/all_mod_models_I_only", header=FALSE)
# #'split.time', 'theta', 'rho', 'log.likelihood'
# 
# ggplot(ImchDat_I_only, aes(x=V2)) + geom_histogram(binwidth =0.01) + xlim(0,0.2)    #Theta
# ggplot(ImchDat_I_only, aes(x=V1)) + geom_histogram(binwidth =0.001) + xlim(0,0.02) + ylim(0,600)    #Split-time



# ImchDat_I_M<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/all_mod_models_I_M", header=FALSE)
# # 'isolation.period', 'migration.period','theta', 'rho', 'migration', 'log.likelihood'
# 
# ggplot(ImchDat_I_M, aes(x=V3)) + geom_histogram(binwidth =0.01) + xlim(0,0.2)   #Theta
# ggplot(ImchDat_I_M, aes(x=V2)) + geom_histogram(binwidth =0.001) + xlim(0,0.02) + ylim(0,600)   #Migration-time
# ggplot(ImchDat_I_M, aes(x=V1)) + geom_histogram(binwidth =0.001) + xlim(0,0.02) + ylim(0,600)   #Split-time


#qplot(ImchDat_a[1], geom="histogram", xlim = c(-6000,1000)) + geom_density()#+ coord_trans(x = "log10")




