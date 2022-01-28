#!/use/bin/env Rscript

if (!require("dplyr")) install.packages("dplyr")
if (!require("data.table")) install.packages("data.table")

library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
      stop("Input arguments are syntfile refreplacerfile qryreplacerfile outfile : ", call.=FALSE)
} 

if (!file.exists(args[1]) || !file.exists(args[2]) || !file.exists(args[3])) {
        stop("Check the existence of syntfile, refreplacerfile and qryreplacerfile", call.=FALSE)
}

syntfile = read.table(args[1], header=F, sep="\t", blank.lines.skip=TRUE)
refreplacer = read.table(args[2], header=F, sep="\t", blank.lines.skip=TRUE) 
qryreplacer = read.table(args[3], header=F, sep="\t", blank.lines.skip=TRUE)
#outfile = read.table(args[4], header=F, sep="\t", blank.lines.skip=TRUE)
outfile = args[4]

if (ncol(syntfile)!=12 || ncol(refreplacer)!=2 || ncol(qryreplacer)!=2) {
    stop("Check validity of of syntfile, refreplacerfile and qryreplacerfile", call.=FALSE)
}

output <- syntfile %>%
    left_join(refreplacer, by = c("V1"="V2")) %>%
    left_join(qryreplacer, by = c("V6"="V2")) %>%
    mutate(V1=V1.y, V6=V1.y.y) %>%
    select(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12)
    
write.table(output, file = outfile, sep="\t", quote = F, row.names = F, col.names=F)    

