library(ggplot2)
library(ggrepel)
library(gridExtra)
library(pracma)
library(plotROC)


# tp_bwa<-909 
# fp_bwa<-26 
# fn_bwa<-16286
# tn_bwa<-52700 
# tpr_bwa<-tp_bwa/(tp_bwa+fn_bwa) 
# fpr_bwa<-fp_bwa/(fp_bwa+tn_bwa) 
# 
# tp_lastz<-7 
# fp_lastz<-3 
# fn_lastz<-17263 
# tn_lastz<-52700 
# tpr_lastz<-tp_lastz/(tp_lastz+fn_lastz) 
# fpr_lastz<-fp_lastz/(fp_lastz+tn_lastz) 
# 
# tp_minimap<-916
# fp_minimap<-12 
# fn_minimap<-16254
# tn_minimap<-52700 
# tpr_minimap<-tp_minimap/(tp_minimap+fn_minimap) 
# fpr_minimap<-fp_minimap/(fp_minimap+tn_minimap) 
# 
# 
# plot(0,0,type="n",ylim=c(0,1),xlim=c(0,1)) 
# points(c(0,fpr_bwa,fpr_lastz,fpr_minimap,1),c(0,tpr_bwa,tpr_lastz,tpr_minimap,1)) 
# lines(c(0,fpr_bwa,fpr_lastz,fpr_minimap,1),c(0,tpr_bwa,tpr_lastz,tpr_minimap,1)) 


data1 <- read.csv("/scr/k61san/nowicklab/SV-detection/BLAT/FINAL/input_ROC.txt", sep="\t", header = TRUE)
plot1 <- ggplot(data1, aes(x=data1$FPR, y=data1$TPR, group=data1$Tool,shape=data1$Tool)) + 
  geom_line() + geom_point(aes(colour=data1$Tool))


data2 <- read.csv("/scr/k61san/nowicklab/SV-detection/BLAT/FINAL/input_PR.txt", sep="\t", header = TRUE)
plot2 <- ggplot(data2, aes(x=data2$Recall, y=data2$Precision, color=factor(data2$Tool),shape=factor(data2$Tool))) + 
  labs(title="(A)                                            ") +
  geom_line() + geom_point() + theme_bw() +
#  theme(axis.title.x = element_text(face = 'bold'),
#                                    axis.title.y = element_text(face = 'bold')) +
  ylab("Precision (TP/(TP+FP))") + xlab("Recall (TP/(TP+FN))") + 
  theme(legend.position = "bottom",legend.title=element_blank(), axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=14,face="bold"))


data3 <- read.csv("/scr/k61san/nowicklab/SV-detection/BLAT/FINAL/input_PPV.txt", sep="\t", header = TRUE)
plot3 <- ggplot(data3, aes(x=data3$PPV, y=data3$Sensitivity, label=data3$Tool)) + 
  labs(title="(B)                                            ") +
  geom_point(aes(shape=factor(data3$Tool)), size=3.5, show.legend=F) +
#  geom_point(aes(color=factor(data3$Threshold)), show.legend=F) + 
  geom_path(aes(color=factor(data3$Tool)), show.legend=T) +
#  geom_point(colour = "black", size = 1.5) + 
  scale_linetype_manual("", values=c(1,2,3,4)) +
  scale_shape_manual("", 
                     values=c(8,9,0,3)) +
  scale_color_discrete(name="") + 
  theme_bw() + #geom_text_repel(point.padding = NA) +
  labs(y = "Sensitivity (TP/(TP+FN))", x= "PPV (TP/(TP+FP))") + #, title="Sensitivity and PPV in relation in called Inversions") +
 theme(legend.position = "bottom", legend.title=element_blank(), axis.text=element_text(size=16,face="bold"),
       axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=14,face="bold")) 
  

svg(filename="/homes/biertank/rohit/Downloads/scripts/simulation_SVs/accuracy_curves.svg", 
    height = 8, width = 12, pointsize = 12)  

grid.arrange(grobs = list(plot2, plot3), ncol=2)

dev.off()


plot3_1 <- 
  ggplot(data3[which(data$Threshold == 0.5) %in% c("Tool","TP","FP","FN","TN","Overlap","Threshold","Sensitivity","PPV")], 
         aes(x=data3$PPV, y=data3$Sensitivity, label=data3$Tool)) + 
  labs(title="(B)                                            ") +
  geom_point(aes(shape=factor(data3_1$Tool)), size=3.5, show.legend=F) + 
  geom_path(aes(color=factor(data3_1$Tool)), show.legend=T) +  scale_linetype_manual("", values=c(1,2,3,4)) +
  scale_shape_manual("", values=c(8,9,0,3)) + scale_color_discrete(name="") + 
  theme_bw() + labs(y = "Sensitivity (TP/(TP+FN))", x= "PPV (TP/(TP+FP))") + 
  theme(legend.position = "bottom", legend.title=element_blank(), axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=14,face="bold")) 

data3_1 <- subset(data3, Threshold == 0.25, select = c("Tool","TP","FP","FN","TN","Overlap","Threshold","Sensitivity","PPV"))
ggplot(data3_1, aes(x=data3_1$PPV, y=data3_1$Sensitivity, label=data3_1$Tool)) + 
  labs(title="0.25 <= Overlap <=1") +
  geom_point(aes(shape=factor(data3_1$Tool), color=factor(data3_1$Tool)), size=3.5, show.legend=T) + 
#  geom_path(aes(color=factor(data3_1$Tool)), show.legend=T) 
#  +  scale_linetype_manual("", values=c(1,2,3,4)) +
  scale_shape_manual("", values=c(8,9,0,3)) + scale_color_discrete(name="") + 
  theme_bw() + labs(y = "Sensitivity (TP/(TP+FN))", x= "PPV (TP/(TP+FP))") + 
  theme(legend.position = "bottom", legend.title=element_blank(), axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5)) 


data2_1 <- subset(data2, Threshold >= 0, select = c("Tool","TP","FP","FN","TN","Overlap","Threshold","Recall","Precision"))
ggplot(data2, aes(x=data2$Recall, y=data2$Precision, color=factor(data2$Tool),shape=factor(data2$Tool))) + 
  labs(title="0.01 <= Overlap <=1") + 
  theme_bw() + geom_point(aes(shape=factor(data2$Tool), show.legend=T)) + 
  geom_path() +
  ylab("Precision (TP/(TP+FP))") + xlab("Recall (TP/(TP+FN))") + 
  theme(legend.position = "bottom",legend.title=element_blank(), axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5))

data2 <- read.csv("/scr/k61san/nowicklab/SV-detection/BLAT/FINAL/input_PR.txt", sep="\t", header = TRUE)
plot_final <- 
  ggplot(data2, aes(x=data2$Overlap, y=data2$Precision, color=factor(data2$Tool),shape=factor(data2$Tool))) + 
  labs(title="(A) Precision at different overlaps") +
  geom_line() + geom_point() + theme_bw() +
          ylab("Precision (TP/(TP+FP))") + xlab("Overlap") + ylim(0,1) + xlim(0,1) +
  theme(legend.position = "bottom",legend.title=element_blank(), axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=14,face="bold"),
        plot.title = element_text(size=15,face="bold", hjust = 0.0)) 

data2 <- read.csv("/scr/k61san/nowicklab/SV-detection/BLAT/FINAL/input_PR.txt", sep="\t", header = TRUE)
newdata <- subset(data2, Overlap == 0.01 | Overlap == 0.09 | Overlap == 0.25 | Overlap == 0.49 | Overlap == 0.77 | Overlap == 0.97 , 
                  select = c("Tool","TP","FP","FN","TN","Overlap","Threshold","Recall","Precision"))
plot_TPs <-  
  ggplot(newdata, aes(x=newdata$Overlap, y=log10(newdata$TP), fill=factor(newdata$Tool))) + 
    labs(title="(B) Number of True positives") +
  geom_bar(stat = "identity", position = "dodge") + 
    theme_bw() +
    ylab("True positives (log10-scaled)") + xlab("Overlap") + 
    scale_x_discrete(labels=c("0.1","0.25","0.5","0.75","1")) + xlim(0,1.1) +
    theme(legend.position = "bottom",legend.title=element_blank(), axis.text=element_text(size=16,face="bold"),
          axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=14,face="bold"),
          plot.title = element_text(size=15,face="bold", hjust = 0.0)) 
  
svg(filename="/homes/biertank/rohit/THESIS/TEX/figures/paper/accuracy_curves_2.svg", 
    height = 8, width = 12, pointsize = 12)   
grid.arrange(grobs = list(plot_final, plot_TPs), ncol=2)
dev.off()

#  scale_shape_manual(values=c(19,20,21))+
#  scale_colour_manual(values=c("blue", "red","gray"))