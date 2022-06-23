
library(ggplot2)
library(gridExtra)

setwd("/homes/biertank/rohit/Downloads/scripts/Genome_annotations")

#col.bil = c("darkblue")
#col.vir = c("darkgreen")

repeat_lbil<-read.table("./repeats_bil_freq.txt", header=F)
repeat_bil <- head(repeat_lbil,15)
repeat_lvir<-read.table("./repeats_vir_freq.txt", header=F)
repeat_vir <- head(repeat_lvir,15)


#png("fig_1_repeat_families.png", width = 11, height = 8, units = 'in', res = 1600)
bitmap("fig_1_repeat_families.tiff", height = 8, width = 12, units = 'in', res=1600)
par(mar=c(2,6,4,2)+0.1, mgp=c(4,1,0))
par(oma=c(6,4,0,0) )
barplot(rbind(repeat_bil$V2,repeat_vir$V2), col=c("blue2","green3"),
        xlim=c(100,250000), cex.names=0.8, xaxt = "n", xpd = FALSE,  
        beside = TRUE,names.arg = repeat_vir$V1, horiz=TRUE, las=1, space=c(0,1),
        main="Repeat-family abundance in Lacertids", col.main=c("red4"))
legend("topright", legend = c("L. bilineata", "L. viridis"), cex=1,
       fill = c("blue2", "green3"), pt.cex = 1)
mtext("Type of Repeat-family",side=2,line=7)
mtext("Number of elements",side=1,line=2)
axis(1, cex.axis=0.8, c("2500","50000","90000","200000"), las = 1)
box()
#dev.copy(png,'fig_1_repeat_families.png')
dev.off()


full_gene_lbil<-read.table("./infile_gf_full_lineata_freq.txt", header=F)
fg_bil <- head(full_gene_lbil,15)
full_gene_lvir<-read.table("./infile_gf_full_viridis_freq.txt", header=F)
fg_vir <- head(full_gene_lvir,15)


bitmap("fig_2_annot_gene_families.tiff", height = 8, width = 12, units = 'in', res=1600)
par(mar=c(2,6,4,2)+0.1, mgp=c(4,1,0))
par(oma=c(6,4,0,0) )
barplot(rbind(fg_bil$V2,fg_vir$V2), col=c("blue2","green3"),
        xlim=c(10,750), cex.names=0.8, xpd = FALSE, #xaxt = "n",
        beside = TRUE,names.arg = fg_vir$V1, horiz=TRUE, las=1, space=c(0,1),
        main="Gene-family abundance in Lacertids", col.main=c("indianred4"))
legend("topright", legend = c("L. bilineata", "L. viridis"), cex=1,
       fill = c("blue2", "green3"), pt.cex = 1)
mtext("Type of pFAM domain",side=2,line=5)
mtext("Number of elements",side=1,line=2)
#axis(1, cex.axis=0.8, c("50","100","150","200","500"), las = 1)
box()
#dev.copy(png,'fig_1_repeat_families.png')
dev.off()


micro_lbil<-read.table("infile_micro_rfam_bil_ovlrm.bed", header=F)
mic_bil <- head(micro_lbil,10)
micro_lvir<-read.table("infile_micro_rfam_bil_ovlrm.bed", header=F)
mic_vir <- head(micro_lvir,10)


bitmap("fig_3_microrna.tiff", height = 8, width = 12, units = 'in', res=1600)
par(mar=c(2,6,4,2)+0.1, mgp=c(4,1,0))
par(oma=c(6,4,0,0) )
barplot(rbind(mic_bil$V2,mic_vir$V2), col=c("blue2","green3"),
        xlim=c(10,6000), cex.names=0.8, xpd = FALSE, xaxt = "n",   
        beside = TRUE,names.arg = mic_vir$V1, horiz=TRUE, las=1, space=c(0,1),
        main="microRNA abundance in Lacertids", col.main=c("violetred4"))
legend("topright", legend = c("L. bilineata", "L. viridis"), cex=1,
       fill = c("blue2", "green3"), pt.cex = 1)
mtext("Gene-family of microRNA",side=2,line=5)
mtext("Number of elements",side=1,line=2)
axis(1, cex.axis=0.8, c("30","350","1000","2000","3000","4500"), las = 1)
box()
#dev.copy(png,'fig_1_repeat_families.png')
dev.off()


tran_bil<-read.table("infile_ilineata_transcript_freq.txt", header=F)
tran_vir<-read.table("infile_vir_transcript_freq.txt", header=F)

bitmap("fig_4_transcripts.tiff", height = 8, width = 12, units = 'in', res=1600)
par(mar=c(2,6,4,2)+0.1, mgp=c(4,1,0))
par(oma=c(6,4,0,0) )
barplot(rbind(tran_bil$V2,tran_vir$V2), col=c("midnightblue","springgreen4"),
        xlim=c(1,65000), cex.names=0.8, xpd = FALSE, #xaxt = "n",   
        beside = TRUE,names.arg = tran_vir$V1, horiz=TRUE, las=1, space=c(0,1),
        main="de-novo transcripts assembled in Lacertids", col.main=c("dodgerblue4"))
legend("topright", legend = c("L. bilineata", "L. viridis"), cex=1,
       fill = c("blue2", "green3"), pt.cex = 1)
mtext("Source / Tissue",side=2,line=5)
mtext("Number of transcripts (including isoforms)",side=1,line=2)
#axis(1, cex.axis=0.8, c("100","1000","2000","3000","4500","6000"), las = 1)
box()
#dev.copy(png,'fig_1_repeat_families.png')
dev.off()



sno_lbil<-read.table("infile_micro_rfam_bil_ovlrm.bed", header=F)
sno_bil <- head(micro_lbil,8)
sno_lvir<-read.table("infile_micro_rfam_vir_ovlrm.bed", header=F)
sno_vir <- head(micro_lvir,8)


bitmap("fig_5_snorna.tiff", height = 8, width = 12, units = 'in', res=1600)
par(mar=c(2,6,4,2)+0.1, mgp=c(4,1,0))
par(oma=c(6,4,0,0) )
barplot(rbind(sno_bil$V2,sno_vir$V2), col=c("blue2","green3"),
        xlim=c(10,6000), cex.names=0.8, xpd = FALSE, xaxt = "n",   
        beside = TRUE,names.arg = sno_vir$V1, horiz=TRUE, las=1, space=c(0,1),
        main="microRNA abundance in Lacertids", col.main=c("violetred4"))
legend("topright", legend = c("L. bilineata", "L. viridis"), cex=1,
       fill = c("blue2", "green3"), pt.cex = 1)
mtext("Gene-family of snoRNA",side=2,line=5)
mtext("Number of elements",side=1,line=2)
axis(1, cex.axis=0.8, c("100","1000","2000","3000","4500","6000"), las = 1)
box()
dev.off()


#   geom_bar(position = "dodge", col="black", fill=col.bil, stat="identity",
#            aes(y=repeat_bil$V2)) + coord_flip() +
#   ggtitle("Top-15 Repeat-families in Lacertids") + 
#   labs(x="Type of Repeat-family", y="Number of elements") 
# 
# rep_lv <- ggplot(repeat_vir, aes(x = repeat_vir$V1)) + theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_bar(position = "dodge", col="black", fill=col.vir, stat="identity",
#            aes(y=repeat_vir$V2)) + coord_flip() +
#   ggtitle("Top-15 Repeat-families in Lacertids") + 
#   labs(x="Type of Repeat-family", y="Number of elements") 
# 
# grid.arrange(rep_lb,rep_lb, ncol = 1)
