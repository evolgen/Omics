setwd("~/Downloads/scripts/population_estimates/")
mu <- 0.8e-9
gen <- 1

bilDat_1<-read.table("/scr/k61san/nowicklab/Lacerta/MSMC/input_1/Estimate_heps/bil.msmc.final.txt", header=TRUE)
bilDat_3<-read.table("/scr/k61san/nowicklab/Lacerta/MSMC/input_1/estimate_heps/cross.msmc.final.txt", header=TRUE)

plot(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, log="x",ylim=c(0,10000000),
     type="n", main = "Plotting MSMC for bilineata and bil+adr", xlab="Years ago", ylab="MSMC2 Effective population size")

lines(bilDat_1$left_time_boundary/mu*gen, (1/bilDat_1$lambda)/mu, type="s", col="blue")

lines(bilDat_3$left_time_boundary/mu*gen, (1/bilDat_3$lambda)/mu, type="s", col="red")

legend("topright",legend=c("bilineata","bil+adr"), col=c("blue","red"), lty=c(1,1))



