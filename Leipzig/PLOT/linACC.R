library(rphast)
library(linACC)

ali_file = system.file("extdata", "/scr/k61san/nowicklab/Lacerta/Multiz/Tetrapods/new_mod_mafs/Lvir_1.maf", package="linACC")
ali = read.msa(ali_file)
ali = read.msa("/scr/k61san/nowicklab/Lacerta/Multiz/Tetrapods/new_mod_mafs/Lvir_1.maf")
# ali_file = system.file("extdata", "/scr/k61san/nowicklab/Lacerta/Multiz/Selected.maf", package="linACC")
# ali_fast = read.msa(ali_file)

mod_file = system.file("extdata", "/scr/k61san/nowicklab/Lacerta/Multiz/Tetrapods/Lacerta_init.mod", package="linACC")
mod      = read.tm(mod_file)
mod = read.tm("/scr/k61san/nowicklab/Lacerta/Multiz/Tetrapods/Lacerta_init.mod")

##- what we looking at
#====================
NBITS    = 14
specs    = c("xenTro3","anoCar2","lacvir1","lacbil1","galGal3","allMis1","hg19")
specs    = as.vector(rbind(specs,specs))
specs = paste(specs,c(".sel",".bgc"),sep="")
PCUT     = 0.005 #- stop traversing if all over sig-thresh


#- enumerate the models for five leaf branches
#=============================================
mod.vecs = int2bitVec(0:(2^NBITS-1),NBITS)
mod.strs = int2bitStr(0:(2^NBITS-1),NBITS)
mod.vecs = t(apply(mod.vecs,1,as.logical))

#- reorder by numbers of parameters
ord            = order(apply(mod.vecs,2,sum))
mod.vecs       = mod.vecs[,ord]
mod.strs       = mod.strs[ord]
rownames(mod.vecs) = specs
colnames(mod.vecs) = mod.strs

#- make the adjacency matrix (nested *with one extra par*)
#=========================================================
amat = getAdjacency(mod.vecs)

ftd = fitNested(ali,mod.vecs,mod,amat)
cla = classifyNested(ftd,mod.vecs,pcut=.01)
cla


#- Alignment gets classified as accelerated in gorilla and gibbon
# ftd = fitNested(ali_fast,mod.vecs,mod,amat)
# cla = classifyNested(ftd,mod.vecs,pcut=.01)
# cla


