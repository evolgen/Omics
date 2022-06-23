setwd("~/Downloads/scripts/genome_metrics/")
library(ggplot2)
library(reshape)
library(lattice)
library(gridExtra)
library(bmp)

dat<-read.table("~/Downloads/scripts/genome_metrics/N50_stats_old.txt", sep="\t", header=TRUE)

dfm <- melt(dat[,c('Species','Pure.illumina','EC.illumina','EC.pacbio','Scaffolded.EC.pacbio','Hybrid')],id.vars = 1)

svg("/homes/biertank/rohit/Downloads/scripts//genome_metrics/N50_values_genome_paper.svg", height = 8, width = 12, pointsize=12)

#bitmap("/homes/biertank/rohit/Downloads/scripts/genome_metrics/N50_values_genome.tiff", height = 8, width = 12, units = 'in', res=1600)

ggplot(dfm, aes(x = Species,y = log10(value))) + 
  geom_bar(aes(fill = variable), stat="identity", color="black", position = "dodge") + 
  geom_text(aes(label=value, group=variable), vjust=-1.1, position=position_dodge(width=0.9), size=4) +
  scale_y_continuous() + ylab("Genomic N50-values (bp)") + xlab("Name of the species") + 
  theme_bw() + 
  theme(axis.text=element_text(size=12, face = "bold.italic"), axis.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=13), legend.key.size = unit(1.25,"line")) +
  theme(legend.title=element_blank())


dev.off()


  # theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(face = 'bold'),
  #                                            axis.title.y = element_text('bold'))


spec_names <- c("Pure.illumina"="Illumina-only","EC.illumina"="Error-corrected Illumina",
                "EC.pacbio"="Error-corrected PacBio","Scaffolded.EC.pacbio"="Scaffolded error-corrected PacBio",
                "Hybrid"="Hybrid-assembly")


dat<-read.table("~/Downloads/scripts/genome_metrics/N50_stats.txt", sep="\t", header=TRUE)

dfm <- melt(dat[,c('Species','Illumina','Illumina.Corrected','PacBio.Corrected','PacBio.Corrected.Extended','Hybrid')],id.vars = 1)


svg("/homes/biertank/rohit/Downloads/scripts//genome_metrics/N50_values_genome_nolog.svg", height = 8, width = 12, pointsize=12)

ggplot(dfm, aes(x = Species,y = value)) + 
  geom_bar(aes(fill = variable), stat="identity", color="black", position = "dodge") + 
  geom_text(aes(label=value, group=variable), vjust=-1.1, position=position_dodge(width=0.9), size=5) +
  scale_y_continuous() + ylab("Genomic N50-values (log-scaled)") + xlab("Name of the species") + 
  theme_bw() + 
  theme(axis.text=element_text(size=14, face = "bold.italic"), axis.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=18), legend.key.size = unit(1.25,"line")) +
  theme(legend.title=element_blank())


dev.off()


# theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(face = 'bold'),
#                                            axis.title.y = element_text('bold'))



# dat2<-read.table("~/Downloads/scripts/genome_metrics/N50_stats.txt", sep="\t", header=TRUE)
# dfm2 <- melt(dat2, id.vars = "Species")

# bitmap("N50_values_genome.tiff", height = 8, width = 12, units = 'in', res=1600)
# ggplot(dfm2,aes(x=variable,y=value,fill=factor(Species))) +
#  geom_bar(stat="identity",position="dodge") + 
#  xlab("Name of the species")+ylab("Genomic N50-values (log-scaled)") + labs(title="Contiguity of Genomic assemblies") +
#  theme_bw() + scale_fill_manual(values = c("blue2","green3")) + theme(legend.title=element_blank()) +
#  geom_text(aes(label=value), vjust=-0.5, position=position_dodge(width=1), size=4)

#dev.off()
#par(mfrow = c(1,1))



