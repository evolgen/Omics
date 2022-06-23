library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(reshape)


setwd("/scr/k61san/nowicklab/Lacerta/RepeatMasker/")

dat2 <- read.table(file = "data_repeatmasker_Tetrapoda.txt_filt_2",header = F,sep = "\t")
dat2 <- dat2[,c(1,2,5)]

dfm2 <- melt(dat2, id.vars=c(1,2),measure.vars = c(3))

dfm2$V1 <- factor(dfm2$V1, levels = c("lepOcu1","xenTro3","cheMyd","allMis1",
                                      "pytBit5","anoCar2","gekJap1","lacVir1","lacBil1",
                                      "taeGut1","galGal4","ornAna5","mm9","hg38"))


dfm2 <- dfm2[order(dfm2$V1),]

#spec_order <- c("asp","cmy","psi","cpc","ggn","cpr","ams","asi","pbi","tsi","vbe",
#                "oha","pmu","cmc","chr",acr","lvi","lbi","gja")

spec_names2 <- c("Lepisosteus oculatus","Xenopus tropicalis","Chelonia mydas",
                "Alligator mississippiensis","Python bivittatus",
                "Anolis carolinensis","Lacerta viridis","Lacerta bilineata","Gekko japonicus",
                "Taeniopygia guttata","Gallus gallus","Ornithorhynchus anatinus",
                "Mus musculus", "Homo sapiens")

# rect <- data.frame(xmin="lacVir1", xmax="lacBil1", ymin=-Inf, ymax=Inf)
# 
# bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/Repeatmasker_species_tetrapods_2.tiff", 
#        height = 8, width = 12, units = 'in', res=1200)

# ggplot(dfm2, aes(x=dfm2$V1,y=log10(value))) +
#   geom_line(aes(group=factor(dfm2$V2), linetype = factor(dfm2$V2), 
#                 color = factor(dfm2$V2)), size=1)

# col_pal <- (palette(gray(seq(0,.9,len = 6))))
# colfunc <- colorRampPalette(c("black", "white"))
# col_pal <- colfunc(6)

# start <- dfm2$V1 == "lacVir1"
# end <- dfm2$V1 == "lacBil1"

svg(filename = "/homes/biertank/rohit/Downloads/scripts/Genome_annotations/Repeatmasker_species_tetrapods_paper.svg", 
    height = 8, width = 12, pointsize = 12)

ggplot(dfm2, aes(x=dfm2$V1,y=log10(value))) +
  geom_line(aes(group=factor(dfm2$V2), linetype = factor(dfm2$V2)), size=0.6) + 
    geom_point(aes(color = factor(dfm2$V1, 
                                labels=spec_names2)), size=2) +
#  scale_size_manual(values = c(0.1,1,0.75,0.9,1.5,0.7)) +
  xlab("Genome versions of Species")+ylab("Percentage of Repeats in the Genome (log-scaled)") + 
  labs(linetype="Type of Repeat") +
  #labs(title="Repeats across Tetrapods predicted with Repeatmasker (Tetrapod lineage)", linetype="Type of Repeat") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate(geom = 'text', x = 8.1, y = 0.94, hjust = 0.35, vjust = -1.2, label = 'L. viridis', fontface = 'bold.italic') +
  annotate(geom = 'text', x = 9.3, y = 0.94, hjust = 0.65, vjust = -1.18, label = 'L. bilineata', fontface = 'bold.italic') +
  theme(axis.text.y=element_text(size=12, face = "bold"), axis.text.x=element_text(size=9, face = "bold"),
	axis.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))


dev.off()


 
# annotate("rect", xmin=7.9, xmax=9.1, ymin=-Inf, ymax=Inf, 
#          fill="yellow", colour="brown", alpha=0.01) #+
#annotate("text", x=1.75, y=17000, label="Region A", size=8)
