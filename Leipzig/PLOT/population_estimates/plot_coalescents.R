setwd("~/Downloads/scripts/population_estimates/")
library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)
library(grid)


# 'isolation.period', 'migration.period','theta', 'rho', 'migration', 'log.likelihood'

#ImchDat_ts<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/I_M_times_scaled", header=FALSE)

dat<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/all_mod_models_I_M", header=FALSE)

ggplot(dat) + geom_density(aes(x=dat$V1), fill="red", alpha=0.1) + geom_density(aes(x=dat$V2), fill = "lightblue", alpha=0.4) +
  xlim(-0.001,0.005) +theme_bw()


dat2<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/I_M_times_scaled", header=FALSE)

ggplot(dat3, aes(x=value, fill=factor(key,labels=c("Isolation-time","Migration-time")), 
                 alpha=0.2)) + geom_density() +
  scale_x_continuous(breaks=c(4000,30000,75000,200000,450000), limits = c(-50000,450000)) + theme_bw() +
  guides(fill=guide_legend(title="Parameter-estimates"), alpha=FALSE) + 
  labs(title="Estimation of Split-times with Isolation-with-Migration model of IMCoalHMM") +
  ylab("Density distribution of estimates") + xlab("Time in Years") + theme(legend.position = c(0.85, 0.9))


dat3 <- gather(dat2)
#bitmap("/homes/biertank/rohit/Downloads/scripts/population_estimates/both_coalascences.tiff", height = 8, width = 12, units = 'in', res=1600)

par(mfrow=c(1,2)) 

plot_1 <- ggplot(dat3, aes(x=value, fill=factor(key,labels=c("Isolation-time","Migration-time")), 
                           alpha=0.2)) + geom_density() +
  scale_x_continuous(breaks=c(4000,30000,75000,200000,450000), limits = c(-50000,450000)) + theme_bw() +
  guides(fill=guide_legend(title="Parameter-estimates"), alpha=FALSE) + 
  labs(title="Split-times estimation with IMCoalHMM") +
  ylab("Density distribution of estimates") + xlab("Time in Years") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                            panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.5, "cm"), legend.position = c(0.75, 0.8))


# plot(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, log="x",ylim=c(0,600000), xlim=c(3000,2000000),
#      type="n", main = "Population history analysis with MSMC2", xlab="Time in years", ylab="Effective population size")
# geom_line(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, type="s", col="red")
# legend("topright",legend=c("mu=1e-8"), col=c("red"), lty=c(1,1))

vp.Right <- viewport(height=unit(0.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","bottom"), 
                           y=0.5, x=0.5)

par(mfrow=c(2,2))
plot(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, log="x",ylim=c(0,600000), xlim=c(3000,2000000),
     type="n", main = "Population history analysis with MSMC2", xlab="Time in years", ylab="Effective population size")
lines(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, type="s", col="red")
legend("topright",legend=c("mu=1e-8"), col=c("red"), lty=c(1,1))


print(plot_1, vp=vp.Right)

#grid.arrange(plot_1,plot_2, ncol=2)
#dev.copy(png,"/homes/biertank/rohit/Downloads/scripts/population_estimates/both_coalescenes.jpeg", height = 8, width = 12, units = 'in', res=1600)
dev.copy(jpeg,"/homes/biertank/rohit/Downloads/scripts/population_estimates/both_coalescenes.jpeg", height = 6, width = 10, units = 'in', res=1000)
dev.off()



#m1 + geom_freqpoly(fill=c("light-blue")) + xlim(0,1.5e+12) + ylim(0,3500)
m1 + geom_area(aes(y = ..count..), fill="lightblue", stat = "bin") + xlim(-0.5,1.5e+12) + ylim(0,3500) +
  geom_smooth(method="gam")
m2 + geom_area(aes(y = ..count..), fill="red", stat = "bin") + xlim(-0.5,1.5e+12) + ylim(0,3500) +
  geom_smooth(method="gam")



mu <- 1e-8
gen <- 1
mu2 <- 0.7e-8

bilDat_1<-read.table("./bil.msmc2.final.txt", header=TRUE)

plot(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, log="x",ylim=c(0,600000), xlim=c(3000,2000000),
     type="n", main = "Population history analysis with MSMC2", xlab="Time in years", ylab="Effective population size")
lines(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, type="s", col="red")
legend("topright",legend=c("mu=1e-8"), col=c("red"), lty=c(1,1))


lines(bilDat_1$left_time_boundary/mu2*gen, (1/bilDat_1$lambda)/mu2, type="s", col="blue")
legend("topright",legend=c("mu=0.7e-8"), col=c("blue"), lty=c(1,1))
