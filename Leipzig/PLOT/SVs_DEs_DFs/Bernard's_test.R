source("https://www.r-statistics.com/wp-content/uploads/2012/01/source_https.r.txt") 
# Making sure we can source code from github
source_https("https://raw.github.com/talgalili/R-code-snippets/master/Barnard.R")

test<-matrix(c(2,2,15,2445),nrow=2, byrow="T")
#LVB_Pos_SVs<-matrix(c(41,2934,30,2715),nrow=2, byrow="T")
fisher.test(test, alternative = "greater")
### Positive results

chisq.test(test)
