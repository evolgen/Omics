library(ggplot2)
library(gridExtra)
library(reshape)
require(grid)
library(ggpubr)
library(cowplot)

setwd("/scr/bloodymary/rohit/Lacerta_viridis/Noncoding/tRNA/")
dat <- read.table(file = "plot_tRNA_species_total.txt",header = T,sep = "\t")
dfm <- melt(dat)[,c(1,2,4)]
dfm$Species <- factor(dfm$Species, levels = c("loc","xtr","pbi","acr","lvi","lbi","gja","ams","gga","tgu","oan","hsa"))


spec_names <- c("loc"="Lepisosteus oculatus", "xtr"="Xenopus tropicalis", "pbi"="Python bivittatus", 
                "acr"="Anolis carolinensis",  "lvi"="Lacerta viridis","lbi"="Lacerta bilineata", "gja"="Gekko japonicus", 
                "ams"="Alligator mississippiensis", "gga"="Gallus gallus", "tgu"="Taeniopygia guttata", 
                "oan"="Ornithorhynchus anatinus", "hsa"="Homo sapiens")

dfm<-dfm[order(dfm2$value),]


types <- unique(dfm$tRNA)
type_1 <- c("Functional-tRNA")
type_2 <- c("Pseudo-tRNA")

dfm1 <- subset(dfm,tRNA %in% type_1)
dfm2 <- subset(dfm,tRNA %in% type_2)



grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + 
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}

plot1 <- ggplot(dfm1, aes(x=Species,y=log(value), label = dfm1$Species)) +
  geom_point(aes(color = factor(dfm1$Species, labels=spec_names)), size=3) +
  ylab("Abundance of functional tRNAs (log-scaled)") + #xlab("Species") + 
  theme_bw() + theme(legend.title=element_blank()) + labs(points="type of tRNA") +
  annotate(geom = 'text', x = 4.7, y = 8.8, hjust = 0.35, vjust = -1.2, label = 'L. viridis', fontface = 'italic') +
  annotate(geom = 'text', x = 6, y = 8.8, hjust = 0.65, vjust = -1.18, label = 'L. bilineata', fontface = 'italic') +
  theme(axis.text=element_text(size=10), axis.title.x=element_blank(),
        axis.title.y=element_text(size=10),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))



plot2 <- ggplot(dfm2, aes(x=Species,y=log(value), label = dfm2$Species)) +
  geom_point(aes(color = factor(dfm2$Species, labels=spec_names)), size=3) +
  ylab("Abundance of pseudo-tRNAs (log-scaled)") + #xlab("Species") + 
  theme_bw() + theme(legend.title=element_blank()) + labs(points="type of tRNA") +
  annotate(geom = 'text', x = 4.7, y = 9.7, hjust = 0.35, vjust = -1.2, label = 'L. viridis', fontface = 'italic') +
  annotate(geom = 'text', x = 6, y = 9.7, hjust = 0.65, vjust = -1.18, label = 'L. bilineata', fontface = 'italic') +
  theme(axis.text=element_text(size=10), axis.title.x=element_blank(),
        axis.title.y=element_text(size=10),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))



figure <- grid_arrange_shared_legend(plot1, plot2, ncol=1,  nrow=2, position = "right") 

svg(filename = "/homes/biertank/rohit/Downloads/scripts/paper/tRNA_types_paper.svg", 
    height = 8, width = 12, pointsize = 12)

figure <- grid_arrange_shared_legend(plot1, plot2, ncol=1,  nrow=2, position = "right") 
# 
# annotate_figure(figure, left = text_grob("Abundance of tRNAs (log-scaled)", size=14,
#                                          color = "black", rot = 90, face = "bold"))

dev.off()


