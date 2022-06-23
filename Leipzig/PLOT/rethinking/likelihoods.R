library(rethinking)
#[rohit@bloodymary Selection]$ printf "\nSV\tType\tGenes_overlap_SVs\tGenes_nonoverlap_SVs\tGenes_overlap_SVs_accl\tGenes_nonoverlap_SVs_accl\n"; for file in /scr/bloodymary/rohit/Lacerta_viridis/Selection/cords_positive_genes.bed /scr/bloodymary/rohit/Lacerta_viridis/Selection/cords_non-positive_cds.bed; do printf "SVs\t"$file"\t"$(bedtools intersect -u -a $file -b <(cat /scr/bloodymary/rohit/Lacerta_viridis/Selection/SV_*.bed) | wc -l)"\t"$(bedtools intersect -v -a $file -b <(cat /scr/bloodymary/rohit/Lacerta_viridis/Selection/SV_*.bed) | wc -l)"\t"$(bedtools intersect -u -a $file -b <(cat /scr/bloodymary/rohit/Lacerta_viridis/Selection/SV_*.bed) | bedtools intersect -u -a - -b /scr/bloodymary/rohit/Lacerta_viridis/Selection/Lacerta_phylop_base-scores.bed_negative_merge | wc -l)"\t"$(bedtools intersect -v -a $file -b <(cat /scr/bloodymary/rohit/Lacerta_viridis/Selection/SV_*.bed) | bedtools intersect -u -a - -b /scr/bloodymary/rohit/Lacerta_viridis/Selection/Lacerta_phylop_base-scores.bed_negative_merge | wc -l)"\n"; done
#SV	Type	Genes_overlap_SVs	Genes_nonoverlap_SVs	Genes_overlap_SVs_accl	Genes_nonoverlap_SVs_accl
#SVs	/scr/bloodymary/rohit/Lacerta_viridis/Selection/cords_positive_genes.bed	192	311	191	309
#SVs	/scr/bloodymary/rohit/Lacerta_viridis/Selection/cords_non-positive_cds.bed	5609	5222	5609	5219

y = 192
n = 192+311
x = 311
k = 2

p = y/n
var = n*p*(1-p)
mu = n*p

likelihood = (p^y)*((1-p)^(n-y))

nL = log(likelihood)
LF = -2*(nL)

AIC = LF + 2*(k)

xi <- rbinom(1000000, 20000,0.3817097)
plot(rbinom(1000000, 20000, 0.3817097), type="p", xlab="Sampling", ylab="Probability")

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

values <- estBetaParams(mu,var)





