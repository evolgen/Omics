setwd("/homes/biertank/rohit/Downloads/scripts/SVs_DEs_DFs/")

library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
#library(reshape2)

data1<-read.table("/scr/bloodymary/rohit/Lacerta_viridis/SVs/plot_SV_counts.txt_paper4", header=F, sep = "\t")
#dfm2<-read.table("filterd_SVs_ranges.txt", header=TRUE, sep = "\t")


#orders <- c("Complex"="Complex","DEL"="Deletion","IDP"="Break-point",
#            "INS"="Insertion","INV"="Inversion","INVDUP"="INV+DUP","TRA"="Translocation")


svg("/homes/biertank/rohit/Downloads/scripts/paper/Size_and_count_SVs_paper4.svg", 
    height = 8, width = 12, pointsize=12)

 ggplot(data1, aes(x=data1$V1,y=log10(data1$V2), fill=data1$V1)) + geom_boxplot(alpha=0.8, show.legend=F) +
  geom_bar(data=data1, aes(y=log10(..count..), fill=data1$V1), alpha=0.5, show.legend=F) +
  scale_y_continuous("Number of rearrangements per category (Bars)", 
                     sec.axis = sec_axis(~., name = "Size distribution of rearrangements per category (Whiskers)")) +
  labs(#title = "Size distribution and the number of genomic rearrangements detected between Lacertids", 
       x="Type of genomic rearrangement") +
  theme_bw()  + 
   theme(axis.text=element_text(size=14, face = "bold"), axis.title=element_text(size=19),
                      legend.text=element_text(size=13, face = "italic"), legend.key.size = unit(1.25,"line"))
 # + geom_blank() # + theme(legend.title=element_blank()) 
  #theme(plot.title = element_text(hjust = 0.5))

dev.off()



