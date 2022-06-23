#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
require(Biostrings)
library("plotly")

setwd("/scr/bloodymary/rohit/Lacerta_viridis/SVs/GC_content")

rearr_del = readDNAStringSet("DEL.fa", "fasta")
rearr_dup = readDNAStringSet("DUP.fa", "fasta")
#rearr_ins = readDNAStringSet("INS.fa", "fasta")
rearr_inv = readDNAStringSet("INV.fa", "fasta")
rearr_tra = readDNAStringSet("TRA.fa", "fasta")

rearr_del[1:3]

letterFrequency(rearr_del[[1]], letters="ACGT", OR=0)
# A    C    G    T 
# 9466 9933 9875 8585 

GC_del = data.frame(ID=names(rearr_del), GC=rowSums(alphabetFrequency(rearr_del)[, c(2,3)]/width(rearr_del))*100)


window = 100
letters = c("G", "C")
# compute the GC content in a sliding window (as a fraction) for a sequence no. 364
gc_del = rowSums(letterFrequencyInSlidingView(rearr_del[[11]], window, letters))/window
plot(gc_del, type = 'l')

#gc_ins = rowSums(letterFrequencyInSlidingView(rearr_ins[[11]], window, letters))/window


gc_del = rowSums(letterFrequencyInSlidingView(rearr_del[[11]], window, letters))/window
plot(1:length(gc_del), gc_del, type="n", xlim=c(0,100000), ylim=c(0,1)) + abline(v=2000)
for (i in 1:length(rearr_del)) {
  gc_del = rowSums(letterFrequencyInSlidingView(rearr_del[[i]], window, letters))/window
  lines(lowess(x = 1:length(gc_del), y= gc_del, f = 0.10), col = 12, lwd = 2) + abline(v=length(gc_del) - 2000)
}

gc_dup = rowSums(letterFrequencyInSlidingView(rearr_dup[[11]], window, letters))/window
plot(1:length(gc_dup), gc_dup, type="n", xlim=c(0,100000), ylim=c(0,1)) + abline(v=2000)
for (i in 1:length(rearr_dup)) {
  gc_dup = rowSums(letterFrequencyInSlidingView(rearr_dup[[i]], window, letters))/window
  lines(lowess(x = 1:length(gc_dup), y= gc_dup, f = 0.10), col = 12, lwd = 2) + abline(v=length(gc_dup) - 2000)
}

gc_inv = rowSums(letterFrequencyInSlidingView(rearr_inv[[11]], window, letters))/window
plot(1:length(gc_inv), gc_inv, type="n", xlim=c(0,100000), ylim=c(0,1)) + abline(v=2000)
for (i in 1:length(rearr_inv)) {
  gc_inv = rowSums(letterFrequencyInSlidingView(rearr_inv[[i]], window, letters))/window
  lines(lowess(x = 1:length(gc_inv), y= gc_inv, f = 0.10), col = 12, lwd = 2) + abline(v=length(gc_inv) - 2000)
}

gc_tra = rowSums(letterFrequencyInSlidingView(rearr_tra[[1]], window, letters))/window
plot(1:length(gc_tra), gc_tra, type="n", xlim=c(0,100000), ylim=c(0,1)) + abline(v=2000)
for (i in 1:length(rearr_tra)) {
  gc_tra = rowSums(letterFrequencyInSlidingView(rearr_tra[[i]], window, letters))/window
  lines(lowess(x = 1:length(gc_tra), y= gc_tra, f = 0.10), col = 12, lwd = 2) + abline(v=length(gc_tra) - 2000)
}


plot_del <- 
  plot(1:length(gc_del), gc_del, type="n", xlim=c(0,10000), ylim=c(0,1)) + lines(lowess(x = 1:length(gc_del), y= gc_del, f = 0.10), col = 12, lwd = 2) +
  abline(v=2000)
plot_dup <- plot(1:length(gc_dup), gc_dup, type="n") + lines(lowess(x = 1:length(gc_dup), y= gc_dup, f = 0.10), col = 12, lwd = 2) +
  abline(v=2000)
#plot_ins <- plot(1:length(gc_ins), gc_ins, type="s") + lines(lowess(x = 1:length(gc_ins), y= gc_ins, f = 0.10), col = 12, lwd = 2) +
#   abline(v=2000)
plot_inv <- plot(1:length(gc_inv), gc_inv, type="n") + lines(lowess(x = 1:length(gc_inv), y= gc_inv, f = 0.10), col = 12, lwd = 2) +
  abline(v=2000)
plot_tra <- plot(1:length(gc_tra), gc_tra, type="n") + lines(lowess(x = 1:length(gc_tra), y= gc_tra, f = 0.10), col = 12, lwd = 2) +
  abline(v=2000)

