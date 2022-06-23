#!/opt/bin/R
library(ggplot2)

test <- read.table("/homes/biertank/rohit/Downloads/pacbio/coverage/jelly.out.sizes", nrows = 5)
classes <- sapply(test,class)


chol <- read.table("/homes/biertank/rohit/Downloads/pacbio/coverage/jelly.out.sizes", colClasses = classes)
x <- chol$V2
normalized = (x-min(x))/(max(x)-min(x))

chol2 <- read.table("/homes/biertank/rohit/Downloads/pacbio/coverage/erct.out.sizes", nrows = 6900000, colClasses = classes)
x <- chol2$V2

chol3 <- read.table("/homes/biertank/rohit/Downloads/pacbio/coverage/pbraw.out.sizes", nrows = 1800000, colClasses = classes)
x <- chol3$V2


figure1 <- ggplot(data=chol, aes(chol$V2)) + geom_histogram(aes(y=..density..),breaks=seq(1, 50000, by=500), col="red", fill="blue") + ggtitle("Histogram for Coverage") + xlab("Contig-length") + ylab("Count")
ggsave(filename="coverage_ctg.jpg", plot=figure1)

figure2 <- ggplot(data=chol2, aes(chol2$V2)) + geom_histogram(aes(y=..density..),breaks=seq(1, 50000, by=500), col="red", fill="blue") + ggtitle("Histogram for Coverage") + xlab("Pacbio-corrected") + ylab("Count")
ggsave(filename="coverage_reads_corr.jpg", plot=figure2)

figure3 <- ggplot(data=chol3, aes(chol3$V2)) + geom_histogram(aes(y=..density..),breaks=seq(1, 50000, by=500), col="red", fill="blue") + ggtitle("Histogram for Coverage") + xlab("Pacbio-raw") + ylab("Count")
ggsave(filename="coverage_reads_raw.jpg", plot=figure3)

#chol <- read.table("/homes/biertank/rohit/Downloads/lacerta/jelly.out.sizes")
figure4 <- ggplot(data=chol, aes(normalized)) + geom_density(aes(y=..density..),breaks=seq(1, 50000, by=500), col="red", fill="blue") + labs(title="Histogram for Coverage") + labs(xlab = "Length", ylab = "Count", xlim=c(100, 50000))
ggsave(filename="coverage_ctg2.jpg", plot=figure4)

figure5 <- ggplot(data=chol, aes(normalized)) + geom_histogram(aes(y=..density..), col="red", fill="green") + labs(title="Histogram for Coverage") + labs(xlab = "Length", ylab = "Count", xlim=c(1, 50000))
ggsave(filename="coverage_ctg3.jpg", plot=figure5)


