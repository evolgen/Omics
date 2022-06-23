setwd("/homes/biertank/rohit/Downloads/scripts")

library(pastecs)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(nortest)


data1<-read.table("/scr/bloodymary/rohit/Lacerta_viridis/Selection/Divergence_collinear_inv.txt", header=T, sep = "\t")

stat.desc(data1$dS)
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 2.202000e+03 0.000000e+00 0.000000e+00 9.513530e-06 6.934290e+00 6.934280e+00 2.907891e+02 1.887355e-02 1.320568e-01 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 1.363547e-02 2.673974e-02 4.094092e-01 6.398509e-01 4.845271e+00 

stat.desc(data1$Omega)
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 2.202000e+03 0.000000e+00 0.000000e+00 1.000000e-03 5.000000e+01 4.999900e+01 7.870500e+03 4.791685e-02 3.574251e+00 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 2.701455e-01 5.297668e-01 1.606989e+02 1.267671e+01 3.546675e+00 


stat.desc(data1[data1$type %in% c("Collinear"),][,"dS"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 2.093000e+03 0.000000e+00 0.000000e+00 9.513530e-06 6.934290e+00 6.934280e+00 2.765145e+02 1.839210e-02 1.321140e-01 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 1.402559e-02 2.750557e-02 4.117293e-01 6.416613e-01 4.856878e+00 


stat.desc(data1[data1$type %in% c("Inverted"),][,"dS"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 1.090000e+02 0.000000e+00 0.000000e+00 7.235160e-05 4.729880e+00 4.729808e+00 1.427456e+01 2.626420e-02 1.309592e-01 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 5.812498e-02 1.152138e-01 3.682579e-01 6.068426e-01 4.633828e+00 


var.test(data1[data1$type %in% c("Collinear"),][,"dS"], data1[data1$type %in% c("Inverted"),][,"dS"])
# F = 1.118, num df = 2092, denom df = 108, p-value = 0.457
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.8350028 1.4464750
# sample estimates:
#   ratio of variances 
# 1.118046 


wilcox.test(data1[data1$type %in% c("Collinear"),][,"dS"], data1[data1$type %in% c("Inverted"),][,"dS"], 
            paired=FALSE, alternative = c("less"))
# V = 2425500, p-value < 2.2e-16
# alternative hypothesis: true location is greater than 0



data2<-read.table("/scr/k61san/nowicklab/Lacerta/FST/Stats_collinear_inv.txt", header=T, sep = "\t")

stat.desc(data2[data2$type %in% c("Collinear"),][,"Dxy"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 7.826000e+03 0.000000e+00 1.000000e+00 1.000000e-01 1.000000e+00 9.000000e-01 6.653186e+03 8.766500e-01 8.501388e-01 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 1.429469e-03 2.802142e-03 1.599151e-02 1.264575e-01 1.487493e-01 

stat.desc(data2[data2$type %in% c("Inverted"),][,"Dxy"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 1.531000e+03 1.000000e+00 7.000000e+00 0.000000e+00 1.000000e+00 1.000000e+00 1.247617e+03 8.333000e-01 8.149033e-01 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 3.658287e-03 7.175788e-03 2.048947e-02 1.431414e-01 1.756545e-01 


var.test(data2[data2$type %in% c("Collinear"),][,"Dxy"], data2[data2$type %in% c("Inverted"),][,"Dxy"])
# F = 0.94794, num df = 7825, denom df = 435, p-value = 0.4271
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.8233426 1.0822039
# sample estimates:
#   ratio of variances 
# 0.9479444 


data3<-read.table("/scr/k61san/nowicklab/Lacerta/FST/Stats_collinear_inv.txt.mpileup", header=T, sep = "\t")

stat.desc(data3[data3$type %in% c("Collinear"),][,"Dxy"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 8.654000e+03 2.969000e+03 0.000000e+00 0.000000e+00 2.500000e-01 2.500000e-01 1.232600e+01 2.000000e-04 1.424312e-03 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 3.763999e-05 7.378334e-05 1.226072e-05 3.501531e-03 2.458401e+00 

stat.desc(data3[data3$type %in% c("Inverted"),][,"Dxy"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 2.307000e+03 8.730000e+02 4.000000e+00 0.000000e+00 5.000000e-01 5.000000e-01 7.150500e+00 5.000000e-04 3.099480e-03 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 3.284889e-04 6.441645e-04 2.489367e-04 1.577773e-02 5.090444e+00 

### No difference in variance through F-test, so that U-test can be used
var.test(data3[data3$type %in% c("Collinear"),][,"Dxy"], data3[data3$type %in% c("Inverted"),][,"Dxy"])
# F test to compare two variances
# F = 0.049252, num df = 8653, denom df = 2306, p-value < 2.2e-16
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.04612552 0.05252601
# sample estimates:
#   ratio of variances 
# 0.04925235 

wilcox.test(data3[data3$type %in% c("Collinear"),][,"Dxy"], data3[data3$type %in% c("Inverted"),][,"Dxy"], 
            +             paired=FALSE, alternative = c("less"))
# W = 9145400, p-value = 1.157e-10
# alternative hypothesis: true location shift is less than 0


 stat.desc(data3[data3$type %in% c("Collinear"),][,"pi_viridis"])
#      nbr.val     nbr.null       nbr.na          min          max        range
# 8.654000e+03 2.761000e+03 0.000000e+00 0.000000e+00 1.468000e-01 1.468000e-01
#          sum       median         mean      SE.mean CI.mean.0.95          var
# 2.630680e+01 5.000000e-04 3.039843e-03 5.724808e-05 1.122199e-04 2.836212e-05
#      std.dev     coef.var
# 5.325610e-03 1.751936e+00
 stat.desc(data3[data3$type %in% c("Inverted"),][,"pi_viridis"])
#      nbr.val     nbr.null       nbr.na          min          max        range
# 2.311000e+03 7.690000e+02 0.000000e+00 0.000000e+00 1.000000e+00 1.000000e+00
#          sum       median         mean      SE.mean CI.mean.0.95          var
# 1.804770e+01 1.900000e-03 7.809476e-03 5.662331e-04 1.110378e-03 7.409526e-04
#      std.dev     coef.var
# 2.722044e-02 3.485566e+00
 stat.desc(data3[data3$type %in% c("Collinear"),][,"pi_bilineata"])
#      nbr.val     nbr.null       nbr.na          min          max        range
# 8.654000e+03 1.271000e+03 0.000000e+00 0.000000e+00 1.639000e-01 1.639000e-01
#          sum       median         mean      SE.mean CI.mean.0.95          var
# 2.358900e+01 1.100000e-03 2.725792e-03 5.451853e-05 1.068693e-04 2.572202e-05
#      std.dev     coef.var
# 5.071688e-03 1.860630e+00
 stat.desc(data3[data3$type %in% c("Inverted"),][,"pi_bilineata"])
#      nbr.val     nbr.null       nbr.na          min          max        range
# 2.307000e+03 5.620000e+02 4.000000e+00 0.000000e+00 1.000000e+00 1.000000e+00
#          sum       median         mean      SE.mean CI.mean.0.95          var
# 2.012660e+01 2.500000e-03 8.724144e-03 6.518832e-04 1.278338e-03 9.803635e-04
#      std.dev     coef.var
# 3.131076e-02 3.588978e+00
 wilcox.test(data3[data3$type %in% c("Collinear"),][,"pi_bilineata"], data3[data3$type %in% c("Inverted"),][,"pi_bilineata"],  paired=FALSE, alternative = c("less"))
#         Wilcoxon rank sum test with continuity correction
# 
# data:  data3[data3$type %in% c("Collinear"), ][, "pi_bilineata"] and data3[data3$type %in% c("Inverted"), ][, "pi_bilineata"]
# W = 8071200, p-value < 2.2e-16
# alternative hypothesis: true location shift is less than 0

 wilcox.test(data3[data3$type %in% c("Collinear"),][,"pi_viridis"], data3[data3$type %in% c("Inverted"),][,"pi_viridis"],  paired=FALSE, alternative = c("less"))
#         Wilcoxon rank sum test with continuity correction
# 
# data:  data3[data3$type %in% c("Collinear"), ][, "pi_viridis"] and data3[data3$type %in% c("Inverted"), ][, "pi_viridis"]
# W = 8567800, p-value < 2.2e-16
# alternative hypothesis: true location shift is less than 0



data4<-read.table("/scr/k61san/nowicklab/Lacerta/FST/Stats_collinear_inv.txt.nowindow", header=T, sep = "\t")


stat.desc(data4[data4$type %in% c("Collinear"),][,"Dxy"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 7.660000e+02 0.000000e+00 0.000000e+00 5.000000e-01 1.000000e+00 5.000000e-01 6.329258e+02 8.216500e-01 8.262739e-01 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 3.882038e-03 7.620712e-03 1.154379e-02 1.074420e-01 1.300320e-01 

stat.desc(data4[data4$type %in% c("Inverted"),][,"Dxy"])
# nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
# 3.390000e+02 3.000000e+00 2.000000e+00 0.000000e+00 1.000000e+00 1.000000e+00 2.719013e+02 8.148000e-01 8.020687e-01 
# SE.mean CI.mean.0.95          var      std.dev     coef.var 
# 7.995453e-03 1.572712e-02 2.167135e-02 1.472119e-01 1.835403e-01 

var.test(data4[data4$type %in% c("Collinear"),][,"Dxy"], data4[data4$type %in% c("Inverted"),][,"Dxy"])
# F = 0.53268, num df = 765, denom df = 338, p-value = 1.65e-12
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.4429296 0.6364899
# sample estimates:
#   ratio of variances 
# 0.5326752 

wilcox.test(data4[data4$type %in% c("Collinear"),][,"Dxy"], data4[data4$type %in% c("Inverted"),][,"Dxy"], 
                         paired=FALSE, alternative = c("two.sided"))
# W = 139120, p-value = 0.05772
# alternative hypothesis: true location shift is not equal to 0



data5<-read.table("/scr/bloodymary/rohit/Lacerta_viridis/Selection/Divergence_collinear_inv.SG.txt.2", header=T, sep = "\t")

 stat.desc(data5[data5$type %in% c("Collinear"),][,"dS"])
#     nbr.val     nbr.null       nbr.na          min          max        range
#8.150000e+02 0.000000e+00 0.000000e+00 2.220100e-05 2.455030e-01 2.454808e-01
#         sum       median         mean      SE.mean CI.mean.0.95          var
#1.699213e+01 1.624920e-02 2.084923e-02 6.829324e-04 1.340516e-03 3.801133e-04
#     std.dev     coef.var
#1.949650e-02 9.351181e-01
 stat.desc(data5[data5$type %in% c("Inverted"),][,"dS"])
#     nbr.val     nbr.null       nbr.na          min          max        range
#5.700000e+01 0.000000e+00 0.000000e+00 3.515240e-03 6.144390e-02 5.792866e-02
#         sum       median         mean      SE.mean CI.mean.0.95          var
#1.348668e+00 2.142840e-02 2.366084e-02 1.932382e-03 3.871027e-03 2.128438e-04
#     std.dev     coef.var
#1.458917e-02 6.165955e-01
 var.test(data5[data5$type %in% c("Collinear"),][,"dS"], data5[data5$type %in% c("Inverted"),][,"dS"])
#data:  data5[data5$type %in% c("Collinear"), ][, "dS"] and data5[data5$type %in% c("Inverted"), ][, "dS"]
#F = 1.7859, num df = 814, denom df = 56, p-value = 0.007499
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
# 1.174759 2.543170
#sample estimates:
#ratio of variances
#           1.78588
 wilcox.test(data5[data5$type %in% c("Collinear"),][,"dS"], 
             data5[data5$type %in% c("Inverted"),][,"dS"], paired=FALSE, alternative = c("less"))
#data:  data5[data5$type %in% c("Collinear"), ][, "dS"] and data5[data5$type %in% c("Inverted"), ][, "dS"]
#W = 19122, p-value = 0.01278
#alternative hypothesis: true location shift is less than 0


 stat.desc(data5[data5$type %in% c("Collinear"),][,"Divergence"])
#      nbr.val     nbr.null       nbr.na          min          max        range
# 8.150000e+02 0.000000e+00 0.000000e+00 1.947460e+03 2.153540e+07 2.153345e+07
#          sum       median         mean      SE.mean CI.mean.0.95          var
# 1.490537e+09 1.425370e+06 1.828880e+06 5.990639e+04 1.175892e+05 2.924852e+12
#      std.dev     coef.var
# 1.710220e+06 9.351185e-01

   stat.desc(data5[data5$type %in% c("Inverted"),][,"Divergence"])
#      nbr.val     nbr.null       nbr.na          min          max        range
# 5.700000e+01 0.000000e+00 0.000000e+00 3.083540e+05 5.389820e+06 5.081466e+06
#          sum       median         mean      SE.mean CI.mean.0.95          var
# 1.183042e+08 1.879680e+06 2.075512e+06 1.695073e+05 3.395639e+05 1.637765e+12
#      std.dev     coef.var
# 1.279752e+06 6.165956e-01

  var.test(data5[data5$type %in% c("Collinear"),][,"Divergence"], data5[data5$type %in% c("Inverted"),][,"Divergence"])
# data:  data5[data5$type %in% c("Collinear"), ][, "Divergence"] and data5[data5$type %in% c("Inverted"), ][, "Divergence"]
# F = 1.7859, num df = 814, denom df = 56, p-value = 0.007499
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#  1.174760 2.543172
# sample estimates:
# ratio of variances
#           1.785881

 wilcox.test(data5[data5$type %in% c("Collinear"),][,"dS"], data5[data5$type %in% c("Inverted"),][,"dS"], paired=FALSE, alternative = c("less"))
# data:  data5[data5$type %in% c("Collinear"), ][, "dS"] and data5[data5$type %in% c("Inverted"), ][, "dS"]
# W = 19122, p-value = 0.01278
# alternative hypothesis: true location shift is less than 0


 
 data6<-read.table("/scr/k61san/nowicklab/Lacerta/FST/Stats_collinear_inv.txt.nowindow.mpileup", header=T, sep = "\t")
 
 stat.desc(data6[data6$type %in% c("Collinear"),][,"Dxy"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 2.903000e+03 6.550000e+02 0.000000e+00 0.000000e+00 1.380000e-02 1.380000e-02 4.051700e+00 4.000000e-04 1.395694e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 3.628615e-05 7.114922e-05 3.822335e-06 1.955079e-03 1.400794e+00 
 
 stat.desc(data6[data6$type %in% c("Inverted"),][,"Dxy"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 8.170000e+02 2.090000e+02 0.000000e+00 0.000000e+00 5.000000e-01 5.000000e-01 2.874800e+00 1.000000e-03 3.518727e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 6.466887e-04 1.269369e-03 3.416745e-04 1.848444e-02 5.253161e+00 
 
 var.test(data6[data6$type %in% c("Collinear"),][,"Dxy"], data6[data6$type %in% c("Inverted"),][,"Dxy"])
 # data:  data6[data6$type %in% c("Collinear"), ][, "Dxy"] and data6[data6$type %in% c("Inverted"), ][, "Dxy"]
 # F = 0.011187, num df = 2902, denom df = 816, p-value < 2.2e-16
 # alternative hypothesis: true ratio of variances is not equal to 1
 # 95 percent confidence interval:
 #   0.01000566 0.01246510
 # sample estimates:
 #   ratio of variances 
 # 0.01118707 
 
 wilcox.test(data6[data6$type %in% c("Collinear"),][,"Dxy"], 
             data6[data6$type %in% c("Inverted"),][,"Dxy"], paired=FALSE, alternative = c("less"))
 # data:  data6[data6$type %in% c("Collinear"), ][, "Dxy"] and data6[data6$type %in% c("Inverted"), ][, "Dxy"]
 # W = 1021900, p-value = 5.566e-10
 # alternative hypothesis: true location shift is less than 0
 
 stat.desc(data6[data6$type %in% c("Collinear"),][,"Fst"])
 # data:  data6[data6$type %in% c("Collinear"), ][, "Dxy"] and data6[data5$type %in% c("Inverted"), ][, "Dxy"]
 # W = 323040, p-value = 0.2721
 # alternative hypothesis: true location shift is less than 0
 
 stat.desc(data6[data6$type %in% c("Inverted"),][,"Fst"])
 # nbr.val      nbr.null        nbr.na           min           max         range           sum        median 
 # 222.00000000    0.00000000    6.00000000   -2.00000000    0.01380000    2.01380000 -194.05080000   -0.68675000 
 # mean       SE.mean  CI.mean.0.95           var       std.dev      coef.var 
 # -0.87410270    0.04711142    0.09284512    0.49272581    0.70194430   -0.80304557 
 
 var.test(data6[data6$type %in% c("Collinear"),][,"Fst"], data6[data6$type %in% c("Inverted"),][,"Fst"])
 # data:  data6[data6$type %in% c("Collinear"), ][, "Fst"] and data6[data6$type %in% c("Inverted"), ][, "Fst"]
 # F = 1.0218, num df = 2807, denom df = 757, p-value = 0.7191
 # alternative hypothesis: true ratio of variances is not equal to 1
 # 95 percent confidence interval:
 #   0.910334 1.142508
 # sample estimates:
 #   ratio of variances 
 # 1.021754 
 
 wilcox.test(data6[data6$type %in% c("Collinear"),][,"Fst"], 
             data6[data6$type %in% c("Inverted"),][,"Fst"], paired=FALSE, alternative = c("less"))
 # data:  data6[data6$type %in% c("Collinear"), ][, "Fst"] and data6[data6$type %in% c("Inverted"), ][, "Fst"]
 # W = 1153000, p-value = 0.9998
 # alternative hypothesis: true location shift is less than 0

 stat.desc(data6[data6$type %in% c("Collinear"),][,"pi_viridis"])
 stat.desc(data6[data6$type %in% c("Inverted"),][,"pi_viridis"])
 wilcox.test(data6[data6$type %in% c("Collinear"),][,"pi_viridis"], 
             data6[data6$type %in% c("Inverted"),][,"pi_viridis"], paired=FALSE, alternative = c("less"))
 
 stat.desc(data6[data6$type %in% c("Collinear"),][,"pi_viridis"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 2.903000e+03 5.040000e+02 0.000000e+00 0.000000e+00 5.810000e-02 5.810000e-02 8.743600e+00 9.000000e-04 3.011919e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 7.990591e-05 1.566781e-04 1.853552e-05 4.305290e-03 1.429418e+00 
 stat.desc(data6[data6$type %in% c("Inverted"),][,"pi_viridis"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 8.170000e+02 1.700000e+02 0.000000e+00 0.000000e+00 5.000000e-01 5.000000e-01 7.822300e+00 3.000000e-03 9.574419e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 1.039804e-03 2.041006e-03 8.833344e-04 2.972094e-02 3.104203e+00 
 wilcox.test(data6[data6$type %in% c("Collinear"),][,"pi_viridis"], 
               +              data6[data6$type %in% c("Inverted"),][,"pi_viridis"], paired=FALSE, alternative = c("less"))
 # data:  data6[data6$type %in% c("Collinear"), ][, "pi_viridis"] and data6[data6$type %in% c("Inverted"), ][, "pi_viridis"]
 # W = 930220, p-value < 2.2e-16
 # alternative hypothesis: true location shift is less than 0
  
 
 stat.desc(data6[data6$type %in% c("Collinear"),][,"pi_bilineata"])
 stat.desc(data6[data6$type %in% c("Inverted"),][,"pi_bilineata"])
 wilcox.test(data6[data6$type %in% c("Collinear"),][,"pi_bilineata"], 
             data6[data6$type %in% c("Inverted"),][,"pi_bilineata"], paired=FALSE, alternative = c("less"))
 
 stat.desc(data6[data6$type %in% c("Collinear"),][,"pi_bilineata"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 2.903000e+03 1.720000e+02 0.000000e+00 0.000000e+00 1.277000e-01 1.277000e-01 7.841400e+00 1.200000e-03 2.701137e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 8.530649e-05 1.672674e-04 2.112570e-05 4.596270e-03 1.701606e+00 
 stat.desc(data6[data6$type %in% c("Inverted"),][,"pi_bilineata"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 8.170000e+02 1.010000e+02 0.000000e+00 0.000000e+00 6.667000e-01 6.667000e-01 9.403000e+00 3.500000e-03 1.150918e-02 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 1.511855e-03 2.967584e-03 1.867423e-03 4.321368e-02 3.754715e+00 
 wilcox.test(data6[data6$type %in% c("Collinear"),][,"pi_bilineata"], 
               +              data6[data6$type %in% c("Inverted"),][,"pi_bilineata"], paired=FALSE, alternative = c("less"))
 # data:  data6[data6$type %in% c("Collinear"), ][, "pi_bilineata"] and data6[data6$type %in% c("Inverted"), ][, "pi_bilineata"]
 # W = 808230, p-value < 2.2e-16
 # alternative hypothesis: true location shift is less than 0

 
 data7<-read.table("/scr/k61san/nowicklab/Lacerta/FST/Stats_background.txt.nowindow.mpileup", header=T, sep = "\t")
 
 stat.desc(data7[data7$type %in% c("Whole"),][,"Dxy"])
 stat.desc(data7[data7$type %in% c("Whole"),][,"Fst"])
 stat.desc(data7[data7$type %in% c("Whole"),][,"pi_viridis"])
 stat.desc(data7[data7$type %in% c("Whole"),][,"pi_bilineata"])
 
 stat.desc(data7[data7$type %in% c("Whole"),][,"Dxy"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 4.597000e+03 1.650000e+02 1.500000e+01 0.000000e+00 1.000000e+00 1.000000e+00 1.362360e+01 1.500000e-03 2.963585e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 2.768860e-04 5.428296e-04 3.524330e-04 1.877320e-02 6.334625e+00 
 
 stat.desc(data7[data7$type %in% c("Whole"),][,"Fst"])
 # nbr.val      nbr.null        nbr.na           min           max         range           sum        median 
 # 4.591000e+03  1.000000e+00  2.100000e+01 -2.000000e+00  1.000000e+00  3.000000e+00 -3.560487e+03 -7.067000e-01 
 # mean       SE.mean  CI.mean.0.95           var       std.dev      coef.var 
 # -7.755363e-01  6.682931e-03  1.310176e-02  2.050413e-01  4.528148e-01 -5.838731e-01 
 
 stat.desc(data7[data7$type %in% c("Whole"),][,"pi_viridis"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 4.598000e+03 5.900000e+01 1.400000e+01 0.000000e+00 3.413000e-01 3.413000e-01 3.135290e+01 3.500000e-03 6.818813e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 2.251947e-04 4.414897e-04 2.331767e-04 1.527012e-02 2.239411e+00 
 
 stat.desc(data7[data7$type %in% c("Whole"),][,"pi_bilineata"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 4.612000e+03 3.100000e+01 0.000000e+00 0.000000e+00 5.333000e-01 5.333000e-01 4.151320e+01 2.900000e-03 9.001127e-03 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 3.622966e-04 7.102748e-04 6.053658e-04 2.460418e-02 2.733456e+00 
 
 
 
 data8<-read.table("/scr/bloodymary/rohit/Lacerta_viridis/Selection/Divergence_Whole.SG.txt.2", header=T, sep = "\t")
 
 stat.desc(data8[data8$type %in% c("Whole"),][,"dS"])
 stat.desc(data8[data8$type %in% c("Whole"),][,"Divergence"])
 
 stat.desc(data8[data8$type %in% c("Whole"),][,"dS"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 2.987000e+03 0.000000e+00 0.000000e+00 1.310330e-05 4.268990e+00 4.268977e+00 7.571341e+01 1.644960e-02 2.534764e-02 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 2.155479e-03 4.226375e-03 1.387788e-02 1.178044e-01 4.647549e+00 

 stat.desc(data8[data8$type %in% c("Whole"),][,"Divergence"])
 # nbr.val     nbr.null       nbr.na          min          max        range          sum       median         mean 
 # 2.987000e+03 0.000000e+00 0.000000e+00 1.149410e+03 3.744730e+08 3.744719e+08 6.641527e+09 1.442950e+06 2.223477e+06 
 # SE.mean CI.mean.0.95          var      std.dev     coef.var 
 # 1.890772e+05 3.707347e+05 1.067858e+14 1.033372e+07 4.647549e+00 

 
 
 data_dense_whole<-read.table("/scr/k61san/nowicklab/Lacerta/FST/Statistics_whole_nowindow.TXT.extend", header=T, sep = "\t")
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"dxy_viridis_bilineata"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"dxy_viridis_bilineata_site"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"Fst_viridis_bilineata"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"Fst_viridis_bilineata_site"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"theta"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"pi_viridis"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"pi_viridis_site"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"pi_bilineata"])
 
 stat.desc(data_dense_whole[data_dense_whole$type %in% c("Whole"),][,"pi_bilineata_site"])
 
 
 
 data_dense_invcol<-read.table("/scr/k61san/nowicklab/Lacerta/FST/Statistics_inv_col_nowindow.TXT.extend", header=T, sep = "\t")
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"dxy_viridis_bilineata"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"dxy_viridis_bilineata_site"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"Fst_viridis_bilineata"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"Fst_viridis_bilineata_site"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"theta"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"pi_viridis"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"pi_viridis_site"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"pi_bilineata"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"pi_bilineata_site"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"dxy_viridis_bilineata"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"dxy_viridis_bilineata_site"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"Fst_viridis_bilineata"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"Fst_viridis_bilineata_site"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"theta"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"pi_viridis"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"pi_viridis_site"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"pi_bilineata"])
 
 stat.desc(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"pi_bilineata_site"])
 
 wilcox.test(data_dense_invcol[data_dense_invcol$type %in% c("Collinear"),][,"theta"], 
             +              data_dense_invcol[data_dense_invcol$type %in% c("Inverted"),][,"theta"], paired=FALSE, alternative = c("less"))
 
 