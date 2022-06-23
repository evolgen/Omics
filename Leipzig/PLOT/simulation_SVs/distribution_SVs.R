setwd("/scr/k61san/nowicklab/SV-detection/BLAT/compare_sim/runs/Manuscript")

library(ggplot2)



mu1 <- 12.5e-9
mu2 <- 6e-9
mu3 <- 2e-9
mu4 <- 1e-9
mu5 <- 5e-10

gen <- 1

bilDat_1<-read.table("./bil.msmc2.final.txt", header=TRUE)

plot(bilDat_1$left_time_boundary/mu4*gen, (1/bilDat_1$lambda)/mu4, log="x",ylim=c(0,5000000), xlim=c(20000,10000000),
     type="n", main = "MSMC2 for population history analysis", xlab="Years ago", ylab="MSMC2 Effective population size")

lines(bilDat_1$left_time_boundary/mu1*gen, (1/bilDat_1$lambda)/mu1, type="s", col="red")
lines(bilDat_1$left_time_boundary/mu2*gen, (1/bilDat_1$lambda)/mu2, type="s", col="blue")
lines(bilDat_1$left_time_boundary/mu3*gen, (1/bilDat_1$lambda)/mu3, type="s", col="orange")
lines(bilDat_1$left_time_boundary/mu4*gen, (1/bilDat_1$lambda)/mu4, type="s", col="purple")
lines(bilDat_1$left_time_boundary/mu5*gen, (1/bilDat_1$lambda)/mu5, type="s", col="black")


legend("topright",legend=c("mu=12.5e-9","mu=6e-9","mu=2e-9","mu=1e-9","mu=5e-9"), col=c("red","blue","orange","purple","black"), lty=c(1,1))




#options(scipen = 999)
ImchDat_a<-read.table("/scr/k61san/nowicklab/Lacerta/IMcoal/all_mod_models_I_M", colClasses="numeric", header=FALSE)

mu=2e-9

ggplot(ImchDat_a) +  geom_density(aes(ImchDat_a$V1),alpha = 0.1,col="red",fill="red")  +
  geom_density(aes(ImchDat_a$V2),alpha = 0.1,col="blue",fill="blue")  +
  xlim(-0.001, 1e-2) + ylab("Split times from syntenic blocks") + xlab("Occurance densities over generations") + theme_bw() +
  title(main = "Estimation of Split-times with IMCoalHMM")

ggplot(ImchDat_a) +  geom_density(aes(ImchDat_a$V1),alpha = 0.1,col="red",fill="red")  +
  geom_density(aes(ImchDat_a$V2),alpha = 0.1,col="blue",fill="blue")  +
  xlim(-0.001, 1e-2) + ylab("Number of Blocks") + xlab("Period of Gene Flow") + theme_bw() +
  scale_fill_manual( values = c("red","blue"),guide = guide_legend(reverse = TRUE))
