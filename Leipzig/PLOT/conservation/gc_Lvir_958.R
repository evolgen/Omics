library("seqinr")
library(gridExtra)
require(grid)

setwd("/scr/k61san/nowicklab/Lacerta/Multiz/Indels/large_Lvir_958")

seq_lacvir1 <- read.fasta(file = "lacvir1.fa")
seq_lacbil1 <- read.fasta(file = "lacbil1.fa")
seq_anoCar2 <- read.fasta(file = "anoCar2.fa2")
seq_galGal3 <- read.fasta(file = "galGal3.fa2")
seq_hg19 <- read.fasta(file = "hg19.fa2")
seq_xenTro3 <- read.fasta(file = "xenTro3.fa2")
seq_allMis1 <- read.fasta(file = "allMis1.fa2")

lvseq <- seq_lacvir1[[1]]
lbseq <- seq_lacbil1[[1]]
acseq <- seq_anoCar2[[1]]
amseq <- seq_allMis1[[1]]
ggseq <- seq_galGal3[[1]]
hgseq <- seq_hg19[[1]]
xtseq <- seq_xenTro3[[1]]

# Then write a function to make a sliding window plot:
  
slidingwindowplot <- function(windowsize, inputseq)
  {
    starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
    n <- length(starts)
    chunkGCs <- numeric(n)
    for (i in 1:n) {
      chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
      chunkGC <- GC(chunk)
      chunkGCs[i] <- chunkGC
    }
    plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}

lv <- slidingwindowplot(50000, lvseq)
lb <- slidingwindowplot(10000, lbseq)
#ac <- slidingwindowplot(2000000, acseq)
gg <- slidingwindowplot(300000, ggseq)
#hg <- slidingwindowplot(1000000, hgseq)
am <- slidingwindowplot(10000, amseq)
xt <- slidingwindowplot(30000, xtseq)
  
grid.arrange(lv, lb, ac, gg, hg, am, ncol=2)

/homes/biertank/rohit/Downloads/scripts/conservation/