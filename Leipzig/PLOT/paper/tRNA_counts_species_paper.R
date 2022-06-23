library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(reshape)
library(plotmath)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/Noncoding/tRNA/")

dat <- read.table(file = "plot_tRNA_species_total.txt",header = T,sep = "\t")

dfm <- melt(dat)[,c(1,2,4)]

dfm$Species <- factor(dfm$Species, levels = c("loc","xtr","pbi","acr","lvi","lbi","gja","ams","gga","tgu","oan","hsa"))


#spec_order <- c("asp","cmy","psi","cpc","ggn","cpr","ams","asi","pbi","tsi","vbe",
#                "oha","pmu","cmc","chr",acr","lvi","lbi","gja")

# spec_names <- c("Apalone spinifera","Chelonia mydas","Pelodiscus sinensis","Chrysemys picta",
#          "Gavialis gangeticus","Crocodylus porosus","Alligator mississippiensis",
#          "Alligator sinensis","Python bivittatus","Thamnophis sirtalis","Vipera berus"
#          ,"Ophiophagus hannah","Protobothrops mucrosquamatus","Crotalus mitchellii",
#          "Crotalus horridus", "Pantherophis guttatus",
#          "Anolis carolinensis","Lacerta viridis","Lacerta bilineata","Gekko japonicus")


spec_names <- c("loc"="Lepisosteus oculatus", "xtr"="Xenopus tropicalis", "pbi"="Python bivittatus", 
                "acr"="Anolis carolinensis",  "lvi"="Lacerta viridis","lbi"="Lacerta bilineata", "gja"="Gekko japonicus", 
                "ams"="Alligator mississippiensis", "gga"="Gallus gallus", "tgu"="Taeniopygia guttata", 
                "oan"="Ornithorhynchus anatinus", "hsa"="Homo sapiens")

dfm<-dfm[order(dfm$value),]


#bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/tRNA_species_tetrapods_paper.tiff", 
#       height = 8, width = 12, units = 'in', res=1200)


svg(filename = "/homes/biertank/rohit/Downloads/scripts/Genome_annotations/tRNA_species_tetrapods_paper.svg", 
    height = 8, width = 12, pointsize = 12)

ggplot(dfm, aes(x=Species,y=log(value), label = dfm$Species)) +
  geom_line(aes(group=tRNA, linetype = factor(tRNA)), size=0.9) + 
  geom_point(aes(color = factor(Species, 
                                labels=spec_names)), size=3) +
  ylab("Total number of tRNAs present (log-scaled)") + # xlab("Species names (abbrv.)")+
  labs(linetype="Type of tRNA") + #title="tRNA counts across Tetrapods and Spotted gar", ) +
  theme_bw() + theme(legend.title=element_blank()) + labs(points="type of tRNA") +
  annotate(geom = 'text', x = 5, y = 9.7, hjust = 0.35, vjust = -1.2, label = 'L. viridis', fontface = 'italic') +
  annotate(geom = 'text', x = 6, y = 9.7, hjust = 0.65, vjust = -1.18, label = 'L. bilineata', fontface = 'italic') +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))
  

dev.off()


# geom_text(data=subset(dfm, dfm$Species == "lvi" & dfm$tRNA == "Pseudo-tRNA"),
#           aes(label="L. viridis"),hjust=0.6,vjust=-1.25) +
# geom_text(data=subset(dfm, dfm$Species == "lbi" & dfm$tRNA == "Pseudo-tRNA"),
#           aes(label="L. bilineata"),hjust=0.75,vjust=-1.25) +

#geom_text(aes(label=ifelse((dfm$Species == "lvi" && dfm$tRNA == "Pseudo-tRNA"), 
#as.character("L. viridis"),'')),hjust=0,vjust=-1.5) +
# annotate("rect", xmin=5.9, xmax=7.1, ymin=-Inf, ymax=Inf, 
#          fill="yellow", colour="sienna2", alpha=0.01) +
#theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
#     axis.title.y = element_text('bold'))


