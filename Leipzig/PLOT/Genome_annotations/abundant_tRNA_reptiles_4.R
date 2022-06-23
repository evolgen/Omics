library(ggplot2)
library(gridExtra)
library(reshape)
require(grid)
library(ggpubr)


setwd("/scr/bloodymary/rohit/Lacerta_viridis/Noncoding/tRNA/")

dat2 <- read.table(file = "abundant_lacertids_allcount.txt2",header = T,sep = "\t")

dfm2 <- melt(dat2)[,c(1,2,4)]

dfm2 <- subset(dfm2,Species %in% c("lvi" , "lbi", "hsa", "gja", "acr", "ams", "cmy", "pbi"))

dfm2$Species <- factor(dfm2$Species, levels = c("cmy", "ams", "pbi", "acr", "gja", "lvi", "lbi", "hsa"))
dfm2 <- dfm2[order(dfm2$Species), ]


spec_order <- c("cmy", "ams", "pbi", "acr", "gja", "lvi", "lbi", "hsa")
dfm2$Species <- ordered(dfm2$Species, spec_order)
dfm2 <- dfm2[with(dfm2, order(tRNA, Species)),]


types <- unique(dfm2$tRNA)
type_1 <- types[1:7]
type_2 <- types[8:14]

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


trnaplot1 <- ggplot(subset(dfm2,tRNA %in% type_1), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity",position="dodge", colour="black") + 
  theme_bw() + theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"),
                    labels = c("C. mydas", "A. mississippiensis", "P. bivittatus", 
                               "A. carolinensis", "G. japonicus", "L. viridis", "L. bilineata",
                                "H. sapiens")) +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14, face = "bold.italic"),
         legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(3,"line"))
        

trnaplot2 <- ggplot(subset(dfm2,tRNA %in% type_2), 
                    aes(x=tRNA,y=log10(value),fill=factor(Species))) +
  geom_bar(stat="identity", position="dodge", colour="black") + 
  #xlab("Type of the tRNA with anti-codon")+ylab("Abundance of tRNA (log-scaled)") + 
  #labs(title="Highly abundant tRNAs in Lacertids - II") +
  theme_bw() + theme(legend.title=element_blank(), axis.title.y=element_blank(),
                     axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("firebrick3","olivedrab","tan4","seagreen3","slateblue1","blue2","green3","orangered"),
                    labels = c("C. mydas", "A. mississippiensis", "P. bivittatus", 
                               "A. carolinensis", "G. japonicus", "L. viridis", "L. bilineata",
                               "H. sapiens")) +
  xlab("Type of the tRNA with anti-codon") + # + ylab("Abundance of tRNA (log-scaled)") +
  theme(axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(3,"line"))


figure <- grid_arrange_shared_legend(trnaplot1, trnaplot2, ncol=1,  nrow=2, position = "right") 


svg(filename = "/homes/biertank/rohit/Downloads/scripts/Genome_annotations/tRNA_abundance_lacertids_paper.svg", 
    height = 8, width = 12, pointsize = 12)


annotate_figure(figure, left = text_grob("Abundance of tRNA (log-scaled)", size=14,
                                         color = "black", rot = 90, face = "bold"))


dev.off()


                # top = text_grob("Visualizing mpg", color = "red", face = "bold", size = 14),
                # bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                #                    hjust = 1, x = 1, face = "italic", size = 10),
#                right = "I'm done, thanks :-)!",
#                fig.lab = "Figure 1", fig.lab.face = "bold"

  

# grid.arrange(trnaplot1, trnaplot2, ncol=1,
#              left = textGrob("Abundance of tRNA (log-scaled)", rot = 90, vjust = 0.75, 
#             top=textGrob("Highly abundant functional tRNAs in Lacertids (tRNA-SCAN)", gp=gpar(fontsize=16,font=8)))

