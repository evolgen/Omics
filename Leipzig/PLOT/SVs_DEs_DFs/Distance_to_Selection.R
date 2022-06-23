library(ggplot2)
library(gridExtra)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/SVs")

data_usgn_neg <- read.table(file = "dist_unsign_neg",header = F,sep = "\t")
data_usgn_relx <- read.table(file = "dist_unsign_relx",header = F,sep = "\t")
data_usgn_pos <- read.table(file = "dist_unsign_pos",header = F,sep = "\t")

data_sgn_neg <- read.table(file = "dist_sign_neg",header = F,sep = "\t")
data_sgn_relx <- read.table(file = "dist_sign_relx",header = F,sep = "\t")
data_sgn_pos <- read.table(file = "dist_sign_pos",header = F,sep = "\t")


P_1 <- ggplot(data_usgn_neg, aes(x=data_usgn_neg$V2, fill=data_usgn_neg$V1)) + xlim(-10000,30000) + 
  geom_density(alpha=.5) + theme_bw() + theme(legend.title=element_blank())

P_2 <- ggplot(data_usgn_relx, aes(x=data_usgn_relx$V2, fill=data_usgn_relx$V1)) + xlim(-10000,30000) + 
  geom_histogram(binwidth=1000, alpha=.5) + theme_bw() + theme(legend.title=element_blank())
  geom_density(alpha=.5) + theme_bw() + theme(legend.title=element_blank())

P_3 <- ggplot(data_usgn_pos, aes(x=data_usgn_pos$V2, fill=data_usgn_pos$V1)) + xlim(-10000,30000) + 
  geom_density(alpha=.5) + theme_bw() + theme(legend.title=element_blank())

P_4 <- ggplot(data_sgn_neg, aes(x=data_sgn_neg$V2, fill=data_sgn_neg$V1)) + xlim(-10000,30000) + 
  geom_density(alpha=.5) + theme_bw() + theme(legend.title=element_blank())

P_5 <- ggplot(data_sgn_relx, aes(x=data_sgn_relx$V2, fill=data_sgn_relx$V1)) + xlim(-10000,30000) + 
  geom_density(alpha=.5) + theme_bw() + theme(legend.title=element_blank())

P_6 <- ggplot(data_sgn_pos, aes(x=data_sgn_pos$V2, fill=data_sgn_pos$V1)) + xlim(-10000,30000) + 
  geom_density(alpha=.5) + theme_bw() + theme(legend.title=element_blank())


grid.arrange(P_1,P_2,P_3,P_4,P_5,P_6,ncol=2)

  ggplot(data_usgn_neg, aes(x=data_usgn_neg$V2, fill=data_usgn_neg$V1)) + xlim(-10000,30000) + 
  #scale_y_log10() + geom_histogram(binwidth=500, alpha=.3) + 
  geom_density(alpha=.5) + theme_bw() + theme(legend.title=element_blank())
  + #   scale_x_continuous(limits = c(0.35, 0.65)) 
  #+ theme_bw() + theme(legend.title=element_blank())

data_usgn_neg <- read.table(file = "test",header = F,sep = "\t")
ggplot(subset(data_usgn_neg, data_usgn_neg$V2<=30000), aes(x=data_usgn_neg$V2, fill=data_usgn_neg$V1)) + xlim(-5000,30000) +
  geom_density(alpha=.1) 
