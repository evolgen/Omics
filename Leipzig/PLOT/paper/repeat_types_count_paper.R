library(ggplot2)
library(gridExtra)
library(reshape)
require(grid)
library(ggpubr)

setwd("/scr/k61san/nowicklab/Lacerta/RepeatMasker/")

dat2 <- read.table(file = "data_repeatmasker_Tetrapoda.txt_filt_2",header = F,sep = "\t")
dat2 <- dat2[,c(1,2,5)]

dfm <- melt(dat2, id.vars=c(1,2),measure.vars = c(3))

dfm$V1 <- factor(dfm$V1, levels = c("lepOcu1","xenTro3","cheMyd","allMis1",
                                      "pytBit5","anoCar2","gekJap1","lacVir1","lacBil1",
                                      "taeGut1","galGal4","ornAna5","mm9","hg38"))


dfm <- dfm2[order(dfm$V1),]

spec_names2 <- c("Lepisosteus oculatus","Xenopus tropicalis","Chelonia mydas",
                 "Alligator mississippiensis","Python bivittatus",
                 "Anolis carolinensis","Lacerta viridis","Lacerta bilineata","Gekko japonicus",
                 "Taeniopygia guttata","Gallus gallus","Ornithorhynchus anatinus",
                 "Mus musculus", "Homo sapiens")


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


types <- unique(dfm$V2)
type_1 <- c("SINEs")
type_2 <- c("LINEs")
type_3 <- c("DNA transposons")

dfm1 <- subset(dfm,dfm$V2 %in% type_1)
dfm2 <- subset(dfm,dfm$V2 %in% type_2)
dfm3 <- subset(dfm,dfm$V2 %in% type_3)


plot1 <- ggplot(dfm1, aes(x=dfm1$V1,y=log10(value))) +
  geom_point(aes(color = factor(dfm1$V1, 
                                labels=spec_names2)), size=2) +
  xlab("Genome versions")+ylab("SINEs") + 
  labs(linetype="Type of Repeat") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate(geom = 'text', x = 8.1, y = 0.43, hjust = 0.35, vjust = -1.2, label = 'L. viridis', fontface = 'italic') +
  annotate(geom = 'text', x = 9.7, y = 0.43, hjust = 0.65, vjust = -1.18, label = 'L. bilineata', fontface = 'italic') +
  theme(axis.text=element_text(size=10), axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))


plot2 <- ggplot(dfm2, aes(x=dfm2$V1,y=log10(value))) +
  geom_point(aes(color = factor(dfm2$V1, 
                                labels=spec_names2)), size=2) +
  xlab("Genome versions")+ylab("LINEs") + 
  labs(linetype="Type of Repeat") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate(geom = 'text', x = 8.1, y = 0.94, hjust = 0.35, vjust = -1.2, label = 'L. viridis', fontface = 'italic') +
  annotate(geom = 'text', x = 9.7, y = 0.95, hjust = 0.65, vjust = -1.18, label = 'L. bilineata', fontface = 'italic') +
  theme(axis.text=element_text(size=10), axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))


plot3 <- ggplot(dfm3, aes(x=dfm3$V1,y=log10(value))) +
  geom_point(aes(color = factor(dfm3$V1, 
                                labels=spec_names2)), size=2) +
  xlab("Genome versions")+ylab("DNA-transposons") + 
  labs(linetype="Type of Repeat") +
  theme_bw() + theme(legend.title=element_blank()) +
  annotate(geom = 'text', x = 8.1, y = 0.22, hjust = 0.35, vjust = -1.2, label = 'L. viridis', fontface = 'italic') +
  annotate(geom = 'text', x = 9.7, y = 0.22, hjust = 0.65, vjust = -1.18, label = 'L. bilineata', fontface = 'italic') +
  theme(axis.text=element_text(size=10), axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))


figure <- grid_arrange_shared_legend(plot1, plot2, plot3, ncol=1,  nrow=3, position = "right") 


svg(filename = "/homes/biertank/rohit/Downloads/scripts/Genome_annotations/repeat_types_paper.svg", 
    height = 8, width = 12, pointsize = 12)

annotate_figure(figure, left = text_grob("Percentage of Repeats in the Genome (log-scaled)", size=14,
                                         color = "black", rot = 90, face = "bold"))

dev.off()


