library(ggplot2)
library(gridExtra)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/Selection/Large_SVs")

collinear <-read.csv(file = "kaks_All_collinear.txt",header = F,sep = "\t")
inverted <-read.csv(file = "kaks_All_inverted.txt",header = F,sep = "\t")

#svg("/scr/bloodymary/rohit/Lacerta_viridis/Wchromosome/Barchart_KaKs_Z.svg", 
#    height = 8, width = 12, pointsize=12)


plot1 <- ggplot(data=collinear, aes(x=collinear$V1,y=collinear$V2), xmax = 3.5, ymax=0.1) + geom_point() +
  labs(title = "Ka/Ks values", x="Ka values", y="Ks values") +
  theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))


plot2 <- ggplot(data=collinear, aes(x=collinear$V1,y=collinear$V2), xmax = 3.5, ymax=0.1) + geom_point() +
  labs(title = "Ka/Ks values", x="Ka values", y="Ks values") +
  theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))

plot3 <- ggplot(data=collinear, aes(x=collinear$V1,y=collinear$V2), xmax = 3.5, ymax=0.1) + geom_point() +
  labs(title = "Ka/Ks values", x="Ka values", y="Ks values") +
  theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))


plot1_b <- ggplot(data=inverted, aes(x=inverted$V1,y=inverted$V2), xmax = 3.5, ymax=0.1) + geom_point() +
  labs(title = "Ka/Ks values", x="Ka values", y="Ks values") +
  theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))

plot2_b <- ggplot(data=inverted, aes(x=inverted$V1,y=inverted$V2), xmax = 3.5, ymax=0.1) + geom_point() +
  labs(title = "Ka/Ks values", x="Ka values", y="Ks values") +
  theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))

plot3_b <- ggplot(data=inverted, aes(x=inverted$V1,y=inverted$V2), xmax = 3.5, ymax=0.1) + geom_point() +
  labs(title = "Ka/Ks values", x="Ka values", y="Ks values") +
  theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))


grid.arrange(plot1, plot2, plot3, plot1_b, plot2_b, plot3_b, ncol=3)

# plot1 <- ggplot(data=collinear, aes(x=collinear$V1)) + geom_bar(stat="bin") +
#   labs(title = "Ka values", x="values", y="Number of Genes") +
#   theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# plot2 <- ggplot(data=collinear, aes(x=collinear$V2)) + geom_bar(stat="bin") +
#   labs(title = "Ks values", x="CAI-values", y="Number of Genes") +
#   theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
# 
# plot3 <- ggplot(data=collinear, aes(x=collinear$V3)) + geom_bar(stat="bin") +
#   labs(title = "Ka-Ks values", x="values", y="Number of Genes") +
#   theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# plot1_b <- ggplot(data=inverted, aes(x=inverted$V1)) + geom_bar(stat="bin") +
#   labs(title = "Inv Ka values", x="values", y="Number of Genes") +
#   theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(data=inverted, aes(x=inverted$V1,y=inverted$V2)) + geom_point() +
#   labs(title = "Inv Ka values", x="values", y="Number of Genes") +
#   theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
# 
# plot2_b <- ggplot(data=inverted, aes(x=inverted$V2)) + geom_bar(stat="bin") +
#   labs(title = "Inv Ks values", x="values", y="Number of Genes") +
#   theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
# 
# plot3_b <- ggplot(data=inverted, aes(x=inverted$V3)) + geom_bar(stat="bin") +
#   labs(title = "Inv Ka-Ks values", x="values", y="Number of Genes") +
#   theme_bw() + theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
# 
# 

