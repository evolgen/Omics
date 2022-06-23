setwd("~/Downloads/scripts/population_estimates/")

library(ggplot2)
library(gridExtra)


mu <- 1e-8
gen <- 1
mu2 <- 0.7e-8

bilDat_1<-read.table("./bil.msmc2.final.txt", header=TRUE)

plot(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, log="x",ylim=c(0,600000), xlim=c(3000,2000000),
     type="n", main = "Population history analysis with MSMC2", xlab="Time (Years ago)", ylab="Effective population size")
lines(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, type="s", col="red")
legend("topright",legend=c("mu=1e-8"), col=c("red"), lty=c(1,1))


lines(bilDat_1$left_time_boundary/mu2*gen, (1/bilDat_1$lambda)/mu2, type="s", col="blue")
legend("topright",legend=c("mu=0.7e-8"), col=c("blue"), lty=c(1,1))


# plot(bilDat_1$left_time_boundary/mu*gen, 2 * bilDat_1$lambda_01 / (bilDat_1$lambda_00 + bilDat_1$lambda_11),
#      xlim=c(1000,500000),ylim=c(0,1), type="n", xlab="Years ago", ylab="Relative cross-coalescence rate")
# lines(bilDat_1$left_time_boundary/mu*gen, 2 * bilDat_1$lambda_01 / (bilDat_1$lambda_00 + bilDat_1$lambda_11), type="s")
# legend("topright",legend=c("Generation-time=1"), col=c("red"), lty=c(1,1))