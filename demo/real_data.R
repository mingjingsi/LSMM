library(LSMM)

Pvalue <- SCZsubset$pval
Z <- SCZsubset[, c(12:20)]
A <- SCZsubset[, 21:147]

fit <- LSMM(Pvalue, Z, A)
assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "local")

