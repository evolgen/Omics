setwd("C:/Users/rohit/Desktop/GenomePaper/SVs/")

library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(reshape2)
library(colorspace)

dfm1<-read.csv("Ratios_of_Accelerated_SVs.txt", header=TRUE, sep = " ")

dat2 <- dfm1[,c(1,2,4,5,6,7,8,9)]

dfm2 <- melt(dat2, id.vars=c(1),measure.vars = c(2,3,4,5,6,7,8))

#dfm2$V1 <- factor(dfm2$SV.type, levels = c("DEL","INS","INV","DUP", "INVDUP","IDP", "TRA","Complex"))


#dfm2 <- dfm2[order(dfm2$SV.type),]

#spec_order <- c("asp","cmy","psi","cpc","ggn","cpr","ams","asi","pbi","tsi","vbe",
#                "oha","pmu","cmc","chr",acr","lvi","lbi","gja")

# spec_names2 <- c("Deletion", "Insertion", "Inversion", "Duplication", "Inverted-duplication",
#                  "Break-point", "Translocation", "Complex-SV")

#rect <- data.frame(xmin="lacVir1", xmax="lacBil1", ymin=-Inf, ymax=Inf)

#bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/Repeatmasker_species_tetrapods_2.tiff", 
#       height = 8, width = 12, units = 'in', res=1200)

pal <- diverge_hcl(7)


ggplot(dfm2, aes(x=dfm2$SV.type,y=log10(value))) +
  geom_line(aes(group=factor(dfm2$variable), linetype = factor(dfm2$variable), 
                color = factor(dfm2$variable)), size=1.5) + 
  geom_point(aes(shape = factor(dfm2$SV.type)), size=2.5) +
  scale_shape_manual(values=seq(0,15)) + scale_color_brewer(palette = "Dark2") +
  xlab("Type of variant overlapping with regions under accelerated evolution")+
  ylab("Frequency of SVs and Accelerated regions (log-scaled)") + 
  labs(title="Accelerated evolution in Lacertids due to SVs") +
   theme_bw() + theme(legend.title=element_blank()) # +
  # annotate("rect", xmin=7.9, xmax=9.1, ymin=-Inf, ymax=Inf, 
  #          fill="yellow", colour="brown", alpha=0.01) #+
#annotate("text", x=1.75, y=17000, label="Region A", size=8)

#dev.off()

