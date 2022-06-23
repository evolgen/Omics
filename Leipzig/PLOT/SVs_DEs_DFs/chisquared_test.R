library(Barnard)
library(Exact)
library(exact2x2)
library(XNomial)

#library(RVAideMemoire)


test<-matrix(c(42,24,140,1468,978,1801),nrow=3, byrow="F")
#fisher.test(test, simulate.p.value = TRUE, B = 1e2)
chisq.test(test, rescale.p = TRUE, simulate.p.value = TRUE, B=1e6)
fisher.test(test)
uncondExact2x2(42,42+140,1468,1468+1801, parmtype="oddsratio", method="wald-unpooled", alternative="greater")
# DELETIONS
X-squared = 52.383, df = NA, p-value = 1e-06
p-value = 4.991e-12
alternative hypothesis: two.sided
data:  x1/n1=(42/182) and x2/n2= (1468/3269)
proportion 1 = 0.23077, proportion 2 = 0.44907, p-value = 1
alternative hypothesis: true p2(1-p1)/[p1(1-p2)] is greater than 1
95 percent confidence interval:
  NA NA
sample estimates:
  p2(1-p1)/[p1(1-p2)] 
2.717009 


test<-matrix(c(7,50,140,129,1744,1801),nrow=3, byrow="F")
#fisher.test(test, simulate.p.value = TRUE, B = 1e2)
chisq.test(test, rescale.p = TRUE, simulate.p.value = TRUE, B=1e6)
fisher.test(test)
uncondExact2x2(7,7+140,129,129+1801, parmtype="oddsratio", method="wald-unpooled", alternative="greater")
# DUPLICATIONS
X-squared = 37.807, df = NA, p-value = 1e-06
p-value = 2.303e-09
alternative hypothesis: two.sided
data:  x1/n1=(7/147) and x2/n2= (129/1930)
proportion 1 = 0.047619, proportion 2 = 0.066839, p-value = 1
alternative hypothesis: true p2(1-p1)/[p1(1-p2)] is greater than 1
95 percent confidence interval:
  NA NA
sample estimates:
  p2(1-p1)/[p1(1-p2)] 
1.432537 


test<-matrix(c(18,47,140,940,1493,1801),nrow=3, byrow="F")
#fisher.test(test, simulate.p.value = TRUE, B = 1e2)
chisq.test(test, rescale.p = TRUE, simulate.p.value = TRUE, B=1e6)
fisher.test(test)
uncondExact2x2(18,18+140,1468,1468+1801, parmtype="oddsratio", method="wald-unpooled", alternative="greater")
# INSERTIONS
X-squared = 54.56, df = NA, p-value = 1e-06
p-value = 9.081e-13
alternative hypothesis: two.sided
data:  x1/n1=(18/158) and x2/n2= (1468/3269)
proportion 1 = 0.11392, proportion 2 = 0.44907, p-value = 1
alternative hypothesis: true p2(1-p1)/[p1(1-p2)] is greater than 1
95 percent confidence interval:
  NA NA
sample estimates:
  p2(1-p1)/[p1(1-p2)] 
6.339688 


test<-matrix(c(8,52,140,77,1778,1801),nrow=3, byrow="F")
#fisher.test(test, simulate.p.value = TRUE, B = 1e2)
chisq.test(test, rescale.p = TRUE, simulate.p.value = TRUE, B=1e6)
fisher.test(test)
uncondExact2x2(42,42+140,1468,1468+1801, parmtype="oddsratio", method="wald-unpooled", alternative="greater")
# INVERSIONS
X-squared = 39.754, df = NA, p-value = 1e-06
p-value = 7.331e-10
alternative hypothesis: two.sided
data:  x1/n1=(42/182) and x2/n2= (1468/3269)
proportion 1 = 0.23077, proportion 2 = 0.44907, p-value = 1
alternative hypothesis: true p2(1-p1)/[p1(1-p2)] is greater than 1
95 percent confidence interval:
  NA NA
sample estimates:
  p2(1-p1)/[p1(1-p2)] 
2.717009 


test<-matrix(c(57,0,140,1817,0,1801),nrow=2, byrow="F")
#fisher.test(test, simulate.p.value = TRUE, B = 1e2)
fisher.test(test)
chisq.test(test, rescale.p = TRUE, simulate.p.value = TRUE, B=1e6)
uncondExact2x2(57,57+140,1817,1817+1801, parmtype="oddsratio", method="wald-unpooled", alternative="greater")
# ANY GR
X-squared = 1160.7, df = NA, p-value = 1e-06
p-value < 2.2e-16
alternative hypothesis: two.sided
proportion 1 = 0.28934, proportion 2 = 0.50221, p-value = 1
alternative hypothesis: true p2(1-p1)/[p1(1-p2)] is greater than 1
95 percent confidence interval:
  NA NA
sample estimates:
  p2(1-p1)/[p1(1-p2)] 
2.477961 



barnard.test(13,294,173,8359)
test <- matrix(c(8,140,77,1801), 2, 2, byrow=F)
exact.test(test, model="multinomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
uncondExact2x2(42,42+140,1468,1468+1801, parmtype="oddsratio", method="wald-unpooled", alternative="greater")



deletion<-matrix(c(42,140,1468,1801),nrow=2, byrow="F")
duplication<-matrix(c(7,140,129,1801),nrow=2, byrow="F")
insertion<-matrix(c(18,140,940,1801),nrow=2, byrow="F")
inversion<-matrix(c(8,140,77,1801),nrow=2, byrow="F")
GR<-matrix(c(57,140,1817,1801),nrow=2, byrow="F")
exact.test(deletion, model="multinomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
exact.test(duplication, model="multinomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
exact.test(insertion, model="multinomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
exact.test(inversion, model="multinomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
exact.test(GR, model="multinomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)

uncondExact2x2(42,42+140,1468,1468+1801, parmtype="oddsratio", method="wald-unpooled", alternative="greater")


deletion<-matrix(c(70,237,3301,5411),nrow=2, byrow="F")
duplication<-matrix(c(12,295,249,8463),nrow=2, byrow="F")
insertion<-matrix(c(33,274,2111,6601),nrow=2, byrow="F")
inversion<-matrix(c(13,294,173,8539),nrow=2, byrow="F")
GR<-matrix(c(95,212,4138,4574),nrow=2, byrow="F")

exact.test(deletion, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
z-unpooled
data:  70 out of 3371 vs. 237 out of 5648
test statistic = -5.8452, first sample size = 3371, second sample size
= 5648, p-value = 1
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion
-0.0211964

exact.test(deletion, model="binomial", alternative="less", method="Z-unpooled", to.plot=FALSE)
data:  70 out of 3371 vs. 237 out of 5648
test statistic = -5.8452, first sample size = 3371, second sample size
= 5648, p-value = 4.554e-08
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion
-0.0211964


> deletion<-matrix(c(3301,5411,70,237),nrow=2, byrow="F")
> exact.test(deletion, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
z-unpooled
data:  3301 out of 3371 vs. 5411 out of 5648
test statistic = 5.8452, first sample size = 3371, second sample size =
  5648, p-value = 4.554e-08
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion 
0.0211964 


> exact.test(duplication, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
z-unpooled
data:  12 out of 261 vs. 295 out of 8758
test statistic = 0.93799, first sample size = 261, second sample size =
  8758, p-value = 0.3131
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion
0.01229352

> exact.test(duplication, model="binomial", alternative="less", method="Z-unpooled", to.plot=FALSE)
data:  12 out of 261 vs. 295 out of 8758
test statistic = 0.93799, first sample size = 261, second sample size =
  8758, p-value = 0.9259
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion
0.01229352


exact.test(insertion, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
z-unpooled
data:  33 out of 2144 vs. 274 out of 6875
test statistic = -6.8822, first sample size = 2144, second sample size
= 6875, p-value = 1
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion
-0.02446275

> exact.test(insertion, model="binomial", alternative="less", method="Z-unpooled", to.plot=FALSE)
data:  33 out of 2144 vs. 274 out of 6875
test statistic = -6.8822, first sample size = 2144, second sample size
= 6875, p-value = 4.754e-07
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion
-0.02446275


> insertion<-matrix(c(2111,6601,33,274),nrow=2, byrow="F")
> exact.test(insertion, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
z-unpooled
data:  2111 out of 2144 vs. 6601 out of 6875
test statistic = 6.8822, first sample size = 2144, second sample size =
  6875, p-value = 4.754e-07
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion 
0.02446275 




exact.test(inversion, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
z-unpooled
data:  13 out of 186 vs. 294 out of 8833
test statistic = 1.9481, first sample size = 186, second sample size =
  8833, p-value = 0.8107
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion
0.0366082

> exact.test(inversion, model="binomial", alternative="less", method="Z-unpooled", to.plot=FALSE)
data:  13 out of 186 vs. 294 out of 8833
test statistic = 1.9481, first sample size = 186, second sample size =
  8833, p-value = 1
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion
0.0366082


exact.test(GR, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)
z-unpooled
data:  95 out of 4233 vs. 212 out of 4786
test statistic = -5.8346, first sample size = 4233, second sample size
= 4786, p-value = 1
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion
-0.02185315

> exact.test(GR, model="binomial", alternative="less", method="Z-unpooled", to.plot=FALSE)
data:  95 out of 4233 vs. 212 out of 4786
test statistic = -5.8346, first sample size = 4233, second sample size
= 4786, p-value = 3.549e-09
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion
-0.02185315





exact.test(GR, model="binomial", alternative="greater", method="Z-unpooled", to.plot=FALSE)


example1<-matrix(c(6,12,12,5),2,2,dimnames=list(c("Group A","Group B"),c("Event","No Event")))


### PAPER
> deletion<-matrix(c(70,237,3301,5411),nrow=2, byrow="F")
> duplication<-matrix(c(12,295,249,8463),nrow=2, byrow="F")
> insertion<-matrix(c(33,274,2111,6601),nrow=2, byrow="F")
> inversion<-matrix(c(13,294,173,8539),nrow=2, byrow="F")
> GR<-matrix(c(95,212,4138,4574),nrow=2, byrow="F")



> exact.test(deletion, model="binomial", alternative="less", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  70 out of 3371 vs. 237 out of 5648
test statistic = -5.3708, first sample size = 3371, second sample size
= 5648, p-value = 7.301e-08
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion 
-0.0211964 


> exact.test(duplication, model="binomial", alternative="greater", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  12 out of 261 vs. 295 out of 8758
test statistic = 1.0793, first sample size = 261, second sample size =
  8758, p-value = 0.1428
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion 
0.01229352 


> exact.test(insertion, model="binomial", alternative="less", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  33 out of 2144 vs. 274 out of 6875
test statistic = -5.4539, first sample size = 2144, second sample size
= 6875, p-value = 4.174e-07
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion 
-0.02446275 


> exact.test(inversion, model="binomial", alternative="greater", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  13 out of 186 vs. 294 out of 8833
test statistic = 2.7248, first sample size = 186, second sample size =
  8833, p-value = 0.0254
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion 
0.0366082 



# [rohit@bloodymary SVs]$ for type in DEL DUP INS INV; do printf $type"\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(fgrep -w $type ALL_SVs.bed_FILT | bedtools intersect -u -b - -a $file2 | wc -l)"\t"$(fgrep -v -w $type ALL_SVs.bed_FILT | bedtools intersect -u -b - -a $file2 | wc -l)"\n"; done; done
# DEL
# Positive        37      24
# Nonpositive     1402    989
# DUP
# Positive        3       49
# Nonpositive     72      1729
# INS
# Positive        18      42
# Nonpositive     915     1461
# INV
# Positive        4       49
# Nonpositive     41      1742


  deletion<-matrix(c(37,24,1402,989),nrow=2, byrow="F")
  duplication<-matrix(c(3,49,72,1729),nrow=2, byrow="F")
  insertion<-matrix(c(18,42,915,1461),nrow=2, byrow="F")
  inversion<-matrix(c(4,49,41,1742),nrow=2, byrow="F")


exact.test(deletion, model="binomial", alternative="less", method="Z-pooled", to.plot=FALSE)
exact.test(duplication, model="binomial", alternative="greater", method="Z-pooled", to.plot=FALSE)
exact.test(insertion, model="binomial", alternative="less", method="Z-pooled", to.plot=FALSE)
exact.test(inversion, model="binomial", alternative="greater", method="Z-pooled", to.plot=FALSE)


> exact.test(deletion, model="binomial", alternative="less", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  57 out of 2708 vs. 27 out of 1562
test statistic = 0.8529, first sample size = 2708, second sample size =
  1562, p-value = 0.8045
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion 
0.003763213 

> exact.test(duplication, model="binomial", alternative="greater", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  3 out of 78 vs. 81 out of 4192
test statistic = 1.206, first sample size = 78, second sample size =
  4192, p-value = 0.1711
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion 
0.01913902 

> exact.test(insertion, model="binomial", alternative="less", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  20 out of 1434 vs. 64 out of 2836
test statistic = -1.9156, first sample size = 1434, second sample size
= 2836, p-value = 0.03109
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
  difference in proportion 
-0.008619994 

> exact.test(inversion, model="binomial", alternative="greater", method="Z-pooled", to.plot=FALSE)

z-pooled

data:  4 out of 45 vs. 80 out of 4224
test statistic = 3.3607, first sample size = 45, second sample size =
  4224, p-value = 0.01568
alternative hypothesis: true difference in proportion is greater than 0
sample estimates:
  difference in proportion 
0.06994949 



###############################################################
### SV against not that particular category or no SV
# [rohit@bloodymary SVs]$ for type in DEL DUP INS INV; do printf $type"\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(fgrep -w $type ALL_SVs.bed_FILT | bedtools intersect -u -a - -b $file2 | wc -l)"\t"$(fgrep -w $type ALL_SVs.bed_FILT | bedtools intersect -v -b - -a $file2 | wc -l)"\n"; done; done
# DEL
# POS	57	160
# NPOS	2651	2216
# DUP
# POS	3	194
# NPOS	75	3546
# INS
# POS	20	179
# NPOS	1414	2703
# INV
# POS	4	193
# NPOS	42	3577
# [rohit@bloodymary SVs]$ printf "ANY\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(cat ALL_SVs.bed_FILT | bedtools intersect -u -a - -b $file2 | wc -l)"\t"$(cat ALL_SVs.bed_FILT | bedtools intersect -v -b - -a $file2 | wc -l)"\n"; done
# ANY
# POS	84	146
# NPOS	4186	1861

> deletion<-matrix(c(57,160,2651,2216),nrow=2, byrow="F")
> duplication<-matrix(c(3,194,75,3546),nrow=2, byrow="F")
> insertion<-matrix(c(20,179,1414,2703),nrow=2, byrow="F")
> inversion<-matrix(c(4,193,42,3577),nrow=2, byrow="F")
> any <-matrix(c(84,146,4186,1861),nrow=2, byrow="F")


### SV against not that particular category
# [rohit@bloodymary SVs]$ for type in DEL DUP INS INV; do printf $type"\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(fgrep -w $type ALL_SVs.bed_FILT | bedtools intersect -u -a - -b $file2 | wc -l)"\t"$(fgrep -v -w $type ALL_SVs.bed_FILT | bedtools intersect -u -b - -a $file2 | wc -l)"\n"; done; done
# DEL
# POS	57	24
# NPOS	2651	989
# DUP
# POS	3	49
# NPOS	75	1729
# INS
# POS	20	42
# NPOS	1414	1461
# INV
# POS	4	49
# NPOS	42	1742

> deletion<-matrix(c(57,24,2651,989),nrow=2, byrow="F")
> duplication<-matrix(c(3,49,75,1729),nrow=2, byrow="F")
> insertion<-matrix(c(20,42,1414,1461),nrow=2, byrow="F")
> inversion<-matrix(c(4,49,42,1742),nrow=2, byrow="F")


# [rohit@bloodymary SVs]$ for type in DEL DUP INS INV; do printf $type"\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(fgrep -w $type ALL_SVs.bed_FILT | bedtools intersect -u -b - -a $file2 | wc -l)"\t"$(fgrep -v -w $type ALL_SVs.bed_FILT | bedtools intersect -u -b - -a $file2 | wc -l)"\n"; done; done
# DEL
# POS	37	24
# NPOS	1402	989
# DUP
# POS	3	49
# NPOS	72	1729
# INS
# POS	18	42
# NPOS	915	1461
# INV
# POS	4	49
# NPOS	41	1742
# [rohit@bloodymary SVs]$ printf "ANY\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(cat ALL_SVs.bed_FILT | bedtools intersect -u -b - -a $file2 | wc -l)"\t"$(cat ALL_SVs.bed_FILT | bedtools intersect -v -b - -a $file2 | wc -l)"\n"; done
# ANY
# POS	51	146
# NPOS	1757	1861


 deletion<-matrix(c(37,24,1402,989),nrow=2, byrow="F")
 insertion<-matrix(c(3,49,72,1729),nrow=2, byrow="F")
 duplication<-matrix(c(3,49,72,1729),nrow=2, byrow="F")
 insertion<-matrix(c(18,42,915,1461),nrow=2, byrow="F")
 inversion<-matrix(c(4,49,41,1742),nrow=2, byrow="F")
 any<-matrix(c(51,146,1757,1861),nrow=2, byrow="F")


> fisher.exact(any,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  any
p-value = 2.7e-10
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.2659 0.5164
sample estimates:
odds ratio 
0.3700622 

> fisher.exact(deletion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  deletion
p-value = 0.7935
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.6302 1.8595
sample estimates:
  odds ratio 
1.087507 

> fisher.exact(duplication,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  duplication
p-value = 0.4656
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.3797 4.6177
sample estimates:
odds ratio 
1.469859 

> fisher.exact(insertion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  insertion
p-value = 0.2259
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.3888 1.1941
sample estimates:
  odds ratio 
0.6844107 

> fisher.exact(inversion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  inversion
p-value = 0.03846
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
1.0940 9.7734
sample estimates:
odds ratio 
3.464223 


[rohit@k61 SVs]$ printf "ALL\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt ; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(cat /scr/bloodymary/rohit/Lacerta_viridis/SVs/HIGH/FINALISED_SVS.BED | bedtools intersect -u -b - -a $file2 | wc -l)"\t"$(cat /scr/bloodymary/rohit/Lacerta_viridis/SVs/HIGH/FINALISED_SVS.BED | bedtools intersect -v -b - -a $file2 | wc -l)"\n"; done
ALL
POS     32      165
NPOS    1177    2441

[rohit@k61 SVs]$ for type in DEL DUP INS INV TRA; do printf $type"\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt ; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(fgrep -w $type /scr/bloodymary/rohit/Lacerta_viridis/SVs/HIGH/FINALISED_SVS.BED | bedtools intersect -u -b - -a $file2 | wc -l)"\t"$(fgrep -v -w $type /scr/bloodymary/rohit/Lacerta_viridis/SVs/HIGH/FINALISED_SVS.BED | bedtools intersect -u -b - -a $file2 | wc -l)"\n"; done; done
DEL
POS     23      14
NPOS    975     362
DUP
POS     0       32
NPOS    6       1171
INS
POS     10      27
NPOS    292     1033
INV
POS     4       28
NPOS    65      1127
TRA
POS     0       32
NPOS    2       1175

[rohit@k61 SVs]$ for type in DEL DUP INS INV TRA; do printf $type"\n"; for file2 in /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_POS.txt /scr/bloodymary/rohit/Lacerta_viridis/Selection/Genes_NPOS.txt ; do printf $(echo $file2 | sed -e 's/.txt//' -e 's/.*Selection\/Genes_//')"\t"$(fgrep -w $type /scr/bloodymary/rohit/Lacerta_viridis/SVs/HIGH/FINALISED_SVS.BED | bedtools intersect -u -b - -a $file2 | wc -l)"\t"$(fgrep -w $type /scr/bloodymary/rohit/Lacerta_viridis/SVs/HIGH/FINALISED_SVS.BED | bedtools intersect -v -b - -a $file2 | wc -l)"\n"; done; done
DEL
POS     23      174
NPOS    975     2643
DUP     
POS     0       197
NPOS    6       3612
INS     
POS     10      187
NPOS    292     3326
INV
POS     4       193
NPOS    65      3553
TRA
POS     0       197
NPOS    2       3616


deletion<-matrix(c(23,14,975,362),nrow=2, byrow="F")
duplication<-matrix(c(0,32,6,1171),nrow=2, byrow="F")
insertion<-matrix(c(10,27,292,1033),nrow=2, byrow="F")
inversion<-matrix(c(4,28,65,1127),nrow=2, byrow="F")
any<-matrix(c(32,165,177,2441),nrow=2, byrow="F")

> fisher.exact(any,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  any
p-value = 1.324e-05
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  1.7669 4.0381
sample estimates:
  odds ratio 
2.673257 

> fisher.exact(deletion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  deletion
p-value = 0.1889
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.3083 1.2159
sample estimates:
odds ratio 
0.6102043 

Warning message:
In mget(objectNames, envir = ns, inherits = TRUE) :
internal error -3 in R_decompress1
> fisher.exact(insertion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  insertion
p-value = 0.4298
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.5818 2.7581
sample estimates:
  odds ratio 
1.30997 

> fisher.exact(inversion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  inversion
p-value = 0.101
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.7730 7.2613
sample estimates:
odds ratio 
2.474228 


deletion2<-matrix(c(23,174,975,2643),nrow=2, byrow="F")
duplication2<-matrix(c(0,197,6,3612),nrow=2, byrow="F")
insertion2<-matrix(c(10,187,292,3326),nrow=2, byrow="F")
inversion2<-matrix(c(4,193,65,3553),nrow=2, byrow="F")

> fisher.exact(any,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  any
p-value = 1.324e-05
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  1.7669 4.0381
sample estimates:
  odds ratio 
2.673257 

> fisher.exact(deletion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  deletion
p-value = 0.1889
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.3083 1.2159
sample estimates:
odds ratio 
0.6102043 

Warning message:
In mget(objectNames, envir = ns, inherits = TRUE) :
internal error -3 in R_decompress1
> fisher.exact(insertion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  insertion
p-value = 0.4298
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.5818 2.7581
sample estimates:
  odds ratio 
1.30997 

> fisher.exact(inversion,alternative = "two.sided",tsmethod = NULL,conf.int = TRUE,conf.level = 0.95,tol = 0.00001)

Two-sided Fisher's Exact Test (usual method using minimum likelihood)

data:  inversion
p-value = 0.101
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.7730 7.2613
sample estimates:
odds ratio 
2.474228 


> exact.test(any, model="binomial", alternative="two.sided", method="boschloo", to.plot=FALSE, cond.row = FALSE)

boschloo

data:  32 out of 197 vs. 177 out of 2618
test statistic = 1.3239e-05, first sample size = 197, second sample
size = 2618, p-value = 1.179e-05
alternative hypothesis: true difference in proportion is not equal to 0
sample estimates:
difference in proportion 
0.09482769 


> exact.test(deletion, model="binomial", alternative="less", method="boschloo", to.plot=FALSE, cond.row = FALSE)

        boschloo

data:  23 out of 37 vs. 975 out of 1337
test statistic = 0.10574, first sample size = 37, second sample size =
1337, p-value = 0.07992
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
difference in proportion 
               -0.107623 

> exact.test(insertion, model="binomial", alternative="less", method="boschloo", to.plot=FALSE, cond.row = FALSE)

        boschloo

data:  10 out of 37 vs. 292 out of 1325
test statistic = 0.82321, first sample size = 37, second sample size =
1325, p-value = 0.7768
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
difference in proportion 
              0.04989291 

> exact.test(inversion, model="binomial", alternative="less", method="boschloo", to.plot=FALSE, cond.row = FALSE)

        boschloo

data:  4 out of 32 vs. 65 out of 1192
test statistic = 0.96986, first sample size = 32, second sample size =
1192, p-value = 0.955
alternative hypothesis: true difference in proportion is less than 0
sample estimates:
difference in proportion 
               0.0704698 



> exact.test(deletion, model="binomial", alternative="two.sided", method="boschloo", to.plot=FALSE, cond.row = FALSE)

boschloo

data:  23 out of 37 vs. 975 out of 1337
test statistic = 0.18886, first sample size = 37, second sample size =
1337, p-value = 0.1836
alternative hypothesis: true difference in proportion is not equal to 0
sample estimates:
difference in proportion 
-0.107623 



> exact.test(insertion, model="binomial", alternative="two.sided", method="boschloo", to.plot=FALSE, cond.row = FALSE)

boschloo

data:  10 out of 37 vs. 292 out of 1325
test statistic = 0.42979, first sample size = 37, second sample size =
1325, p-value = 0.4113
alternative hypothesis: true difference in proportion is not equal to 0
sample estimates:
difference in proportion 
0.04989291 

> exact.test(inversion, model="binomial", alternative="two.sided", method="boschloo", to.plot=FALSE, cond.row = FALSE)

boschloo

data:  4 out of 32 vs. 65 out of 1192
test statistic = 0.10099, first sample size = 32, second sample size =
1192, p-value = 0.09336
alternative hypothesis: true difference in proportion is not equal to 0
sample estimates:
difference in proportion 
0.0704698 


> q<-c(0.05,1.179e-05,0.92,0.24,0.07)
> fdr = fdrtool(q, statistic="pvalue")
> fdr$qval
[1] 2.873088e-02 1.645229e-05 2.626254e-01 8.501329e-02 3.256037e-02
> fdr$lfdr
[1] 0.04883233 0.04883233 1.00000000 1.00000000 0.25251206
