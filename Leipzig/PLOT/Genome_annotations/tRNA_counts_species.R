library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/Noncoding/tRNA/")

dat <- read.table(file = "tRNA_species_total.txt",header = T,sep = "\t")

dfm <- melt(dat)[,c(1,2,4)]

dfm$Species <- factor(dfm$Species, levels = c("asp","cmy","psi","cpc","ggn","cpr","ams","asi","pbi","tsi","vbe",
                                          "oha","pmu","cmc","chr","pgt","acr","lvi","lbi","gja"))

dfm<-dfm[order(dfm$Species),]

#spec_order <- c("asp","cmy","psi","cpc","ggn","cpr","ams","asi","pbi","tsi","vbe",
#                "oha","pmu","cmc","chr",acr","lvi","lbi","gja")

spec_names <- c("Apalone spinifera","Chelonia mydas","Pelodiscus sinensis","Chrysemys picta",
         "Gavialis gangeticus","Crocodylus porosus","Alligator mississippiensis",
         "Alligator sinensis","Python bivittatus","Thamnophis sirtalis","Vipera berus"
         ,"Ophiophagus hannah","Protobothrops mucrosquamatus","Crotalus mitchellii",
         "Crotalus horridus", "Pantherophis guttatus",
         "Anolis carolinensis","Lacerta viridis","Lacerta bilineata","Gekko japonicus")


bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/tRNA_species_reptiles_2.tiff", 
       height = 8, width = 12, units = 'in', res=1200)

ggplot(dfm, aes(x=Species,y=log(value))) +
  geom_line(aes(group=tRNA, linetype = factor(tRNA)), size=0.9) + 
  geom_point(aes(color = factor(Species, 
                                labels=spec_names)), size=3) +
  xlab("Species names (abbrv.)")+ylab("Total number of tRNAs present (log-scaled)") + 
  labs(title="tRNA counts across Reptiles", linetype="Type of tRNA") +
  theme_bw() + theme(legend.title=element_blank()) + labs(points="type of tRNA") +
  annotate("rect", xmin=17.9, xmax=19.1, ymin=-Inf, ymax=Inf, 
           fill="yellow", colour="sienna2", alpha=0.01) #+

dev.off()


#bitmap("/homes/biertank/rohit/Downloads/scripts/Genome_annotations/tRNA_species_reptiles.tiff", height = 20, width = 18, units = 'in', res=2000)

  # ggplot(dfm, aes(x=Species,y=log(value),fill=factor(Species))) +
  # geom_line(aes(group=tRNA, linetype = factor(tRNA)), size=0.9) + 
  # geom_point(aes(color = factor(Species, 
  # labels=c("Alligator mississippiensis","Alligator sinensis","Apalone spinifera",
  #          "Chelonia mydas","Chrysemys picta","Crocodylus porosus","Crotalus mitchellii",
  #          "Gavialis gangeticus","Gekko japonicus","Haliaeetus albicilla","Lacerta viridis",
  #          "Macropus eugenii","Melopsittacus undulatus","Mus musculus","Ophiophagus hannah",
  #          "Pantherophis guttatus", "Protobothrops mucrosquamatus","Python bivittatus",
  #          "Thamnophis sirtalis","Vipera berus"))), size=3) +
  # xlab("Species names (abbrv.)")+ylab("Total number of tRNAs present (log-scaled)") + 
  # labs(title="tRNA counts across Species", linetype="Type of tRNA") +
  # theme_bw() + theme(legend.title=element_blank()) + labs(points="type of tRNA")
    #+ geom_text(aes(label=""), show_guide = FALSE)




  
  grid.newpage()
  
  footnote <- "ams = Alligator mississippiensis\tasi = Alligator sinensis\tacr = Anolis carolinensis\tasp = Apalone spinifera\tcmy = Chelonia mydas\tcpc = Chrysemys picta\tcpr = Crocodylus porosus\tchr = Crotalus horridus\tcmc = Crotalus mitchellii\tgac = Gasterosteus aculeatus\tggn = Gavialis gangeticus\tgja = Gekko japonicus\tlbi = Lacerta bilineata\tlvi = Lacerta viridis\tmdo = Monodelphis domestica\tmmu = Mus musculus\toha = Ophiophagus hannah\tpgt = Pantherophis guttatus\tpsi = Pelodiscus sinensis\tpmu = Protobothrops mucrosquamatus\tpbi = Python bivittatus\ttni = Tetraodon nigroviridis\ttsi = Thamnophis sirtalis\tvbe = Vipera berus\n"
  #footnote <- "acr = Anolis car.\tams = Alligaor miss.\tasi = Alligaor sin.\tasp = Apalone spi.\tchr = Crotalus hor.\tcmc = Crotalus. \tcmy = \tcpc = \tcpr = \tggn = \tgja = \tlbi = \tlvi = \toha = \tpbi = \tpgt = \tpmu = \tpsi = \ttsi = \tvbe = \n"
  
  plot1 <- 
    
cols <- colorRampPalette(brewer.pal(11, "Spectral"))(20)

palette <- rep(c("color1", "color2", ...), length.out = 20)

# or, if you want a random distribution:
# if you want it random, but reproducable,
# backup .Random.seed, or set your own
set.seed(23)
palette <- sample(c("color1", "color2", ...), number, replace = TRUE)

scale_fill_manual(values=palette)


ggplot(dfm, aes(x=Species,y=log(value),fill=factor(Species), show.legend=F)) +
  geom_line(aes(group=tRNA, linetype = factor(tRNA))) + 
  geom_point() +
scale_fill_manual(values=cols,
labels=c("Alligator mississippiensis","Alligator sinensis","Apalone spinifera",
"Chelonia mydas","Chrysemys picta","Crocodylus porosus","Crotalus mitchellii",
"Gavialis gangeticus","Gekko japonicus","Haliaeetus albicilla","Lacerta viridis",
"Macropus eugenii","Melopsittacus undulatus","Mus musculus","Ophiophagus hannah",
"Pantherophis guttatus", "Protobothrops mucrosquamatus","Python bivittatus",
"Thamnophis sirtalis","Vipera berus"))
               
               , color = factor(Species, 
                                               labels=c("Alligator mississippiensis","Alligator sinensis","Apalone spinifera",
                                                        "Chelonia mydas","Chrysemys picta","Crocodylus porosus","Crotalus mitchellii",
                                                        "Gavialis gangeticus","Gekko japonicus","Haliaeetus albicilla","Lacerta viridis",
                                                        "Macropus eugenii","Melopsittacus undulatus","Mus musculus","Ophiophagus hannah",
                                                        "Pantherophis guttatus", "Protobothrops mucrosquamatus","Python bivittatus",
                                                        "Thamnophis sirtalis","Vipera berus"))), size=3) +
  xlab("Species names (abbrv.)")+ylab("Total number of tRNAs present (log-scaled)") + 
  labs(title="tRNA counts in Reptiles and Birds", linetype="Type of tRNA") +
  theme_bw() + theme(legend.title=element_blank(), legend.position='none')
                 
                 
                 
                 ), size=3) +
  scale_fill_manual(values=dfm$Species,
                    labels=c("Alligator mississippiensis","Alligator sinensis","Apalone spinifera",
                                         "Chelonia mydas","Chrysemys picta","Crocodylus porosus","Crotalus mitchellii",
                                         "Gavialis gangeticus","Gekko japonicus","Haliaeetus albicilla","Lacerta viridis",
                                         "Macropus eugenii","Melopsittacus undulatus","Mus musculus","Ophiophagus hannah",
                                         "Pantherophis guttatus", "Protobothrops mucrosquamatus","Python bivittatus",
                                         "Thamnophis sirtalis","Vipera berus")) +
  xlab("Species names (abbrv.)")+ylab("Total number of tRNAs present (log-scaled)") + 
  labs(title="tRNA counts in Reptiles and Birds", linetype="Type of tRNA") +
  theme_bw() + theme(legend.title=element_blank()) #+ #+ guides(color=guide_legend("my title"))
#scale_fill_brewer(palette="Set1")
#+ geom_text(aes(label=""), show_guide = FALSE)




#page1 <- (arrangeGrob(plot1, side = textGrob(footnote, x = 0, hjust = -1, vjust=0.01, gp = gpar(fontface = "italic", fontsize = 12))))

#plot.new()
#page1 <- plot1 + mtext("footnote")
#grid.draw(page1)

#  scale_fill_manual(values = c("firebrick3","olxivedrab","tan4","blue2","green3")) 
#grid.arrange(trnaplot1, trnaplot2, ncol=1)

dev.off()




