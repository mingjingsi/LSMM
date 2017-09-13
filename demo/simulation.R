library(LSMM)
library(pROC)
library(MASS)
library(GPA)

demo(func)
##### 1. Performance of LSMM #####

M      <- 100000   # No. of SNPs
L      <- 10       # No. of fixed effects
K      <- 500      # No. of random effects
Z.perc <- 0.1      # frequency of fixed effects Pr(Z=1)
A.perc <- 0.1      # frequency of random effects Pr(A=1)
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # coefficients of fixed effects
omega  <- 0.2      # proportion of relevant random effects
sigma2 <- 1        # parameter in the spike-slab prior
rep    <- 50       # repeat times

result <- matrix(0, rep, 16)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  data <- generate_data(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2)

  fit <- LSMM(data$Pvalue, data$Z, data$A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")

  result[i, 1:4]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))
  result[i, 5:8]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage1, 1-fit$pi1.stage1))
  result[i, 9:12] <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage2, 1-fit$pi1.stage2))

  relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "global")
  result[i, 13:16]   <- as.numeric(performance(data$eta, relev.Anno, 1-fit$omegak))
}

result1 <- as.data.frame(result)
names(result1) <- c("FDR.SNP", "power.SNP", "AUC.SNP", "pAUC.SNP", "FDR.stage1", "power.stage1",
                    "AUC.stage1", "pAUC.stage1", "FDR.stage2", "power.stage2", "AUC.stage2",
                    "pAUC.stage2", "FDR.Anno", "power.Anno", "AUC.Anno", "pAUC.Anno")

##### 2. Treat all covariates as fixed effects #####

M      <- 100000   # No. of SNPs
L      <- 10       # No. of fixed effects
K      <- 500      # No. of random effects
Z.perc <- 0.1      # frequency of fixed effects Pr(Z=1)
A.perc <- 0.1      # frequency of random effects Pr(A=1)
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # coefficients of fixed effects
omega  <- 0.2      # proportion of relevant random effects
sigma2 <- 1        # parameter in the spike-slab prior
rep    <- 50       # repeat times

result <- matrix(0, rep, 4)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  data <- generate_data(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2)

  fit <- LSMM(data$Pvalue, cbind(data$Z, data$A), NULL)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")

  result[i, 1:4]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))

}

result2 <- as.data.frame(result)
names(result2) <- c("FDR.SNP", "power.SNP", "AUC.SNP", "pAUC.SNP")

##### 3. Adjustment of fixed effects #####

M      <- 100000   # No. of SNPs
L      <- 10       # No. of fixed effects
K      <- 500      # No. of random effects
alpha  <- 0.2      # parameter in the Beta distribution
Z.perc <- 0.1      # frequency of fixed effects Pr(Z=1)
A.perc <- 0.1      # frequency of random effects Pr(A=1)
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # coefficients of fixed effects
omega  <- 0.2      # proportion of relevant random effects
sigma2 <- 1        # parameter in the spike-slab prior
corr   <- 0.2      # correlation coeffecient among fixed effects and the first 50 random effects
rep    <- 50       # repeat times

result <- matrix(0, rep, 8)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  # when fixed effects and random effects are not independent
  data <- generate_data_corr_ZA(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, corr)

  # LSMM
  fit <- LSMM(data$Pvalue, data$Z, data$A)
  relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "global")
  result[i, 1:4]   <- as.numeric(performance(data$eta, relev.Anno, 1-fit$omegak))

  # LSMM without fixed effects
  fit1 <- LSMM(data$Pvalue, NULL, data$A)
  relev.Anno1 <- relev.Anno(fit1, FDRset = 0.1, fdrControl = "global")
  result[i, 5:8]   <- as.numeric(performance(data$eta, relev.Anno1, 1-fit$omegak))
}

result3 <- as.data.frame(result)
names(result3) <- c("FDR.LSMM", "power.LSMM", "AUC.LSMM", "pAUC.LSMM", "FDR.LSRM",
                    "power.LSRM", "AUC.LSRM", "pAUC.LSRM")

##### 4. Simulations based on probit model #####

# 4.1 When the underlying distribution of p-values in non-null group is beta distribution
M      <- 100000   # No. of SNPs
K      <- 500      # No. of random effects
L      <- 10       # No. of fixed effects
Z.perc <- 0.1      # frequency of fixed effects Pr(Z=1)
A.perc <- 0.1      # frequency of random effects Pr(A=1)
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -1       # intercept of the probit model
set.seed(1)
b      <- rnorm(L) # coefficients of fixed effects
omega  <- 0.2      # proportion of relevant random effects
sigma2 <- 1        # parameter in the spike-slab prior
r      <- 5/5      # signal-noise ratio of probit model
rep    <- 50       # repeat times

result <- matrix(0, rep, 16)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  data <- generate_data_probit(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, r)

  fit <- LSMM(data$Pvalue, data$Z, data$A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")

  result[i, 1:4]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))
  result[i, 5:8]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage1, 1-fit$pi1.stage1))
  result[i, 9:12] <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage2, 1-fit$pi1.stage2))

  relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "global")
  out[i, 13:16]   <- as.numeric(performance(data$eta, relev.Anno, 1-fit$omegak))
}

result4.1 <- as.data.frame(result)
names(result4.1) <- c("FDR.SNP", "power.SNP", "AUC.SNP", "pAUC.SNP", "FDR.stage1", "power.stage1",
                    "AUC.stage1", "pAUC.stage1", "FDR.stage2", "power.stage2", "AUC.stage2",
                    "pAUC.stage2", "FDR.Anno", "power.Anno", "AUC.Anno", "pAUC.Anno")

# 4.2 When the underlying distribution of p-values in non-null group is not beta distribution
M      <- 100000   # No. of SNPs
K      <- 500      # No. of random effects
L      <- 10       # No. of fixed effects
Z.perc <- 0.1      # frequency of fixed effects Pr(Z=1)
A.perc <- 0.1      # frequency of random effects Pr(A=1)
alpha  <- 0.2      # parameter in the Beta distribution
beta0  <- -1       # intercept of the probit model
set.seed(1)
b      <- rnorm(L) # coefficients of fixed effects
omega  <- 0.2      # proportion of relevant random effects
sigma2 <- 1        # parameter in the spike-slab prior
r      <- 5/5      # signal-noise ratio of probit model
dist   <- "spiky"  # It can also be "near_normal", "skew" and "big_normal"
rep    <- 50       # repeat times

result <- matrix(0, rep, 16)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  data <- generate_data_probit_dist(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, r, dist)

  fit <- LSMM(data$Pvalue, data$Z, data$A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")

  result[i, 1:4]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))
  result[i, 5:8]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage1, 1-fit$pi1.stage1))
  result[i, 9:12] <- as.numeric(performance(data$gamma, assoc.SNP$gamma.stage2, 1-fit$pi1.stage2))

  relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "global")
  out[i, 13:16]   <- as.numeric(performance(data$eta, relev.Anno, 1-fit$omegak))
}

result4.2 <- as.data.frame(result)
names(result4.2) <- c("FDR.SNP", "power.SNP", "AUC.SNP", "pAUC.SNP", "FDR.stage1", "power.stage1",
                    "AUC.stage1", "pAUC.stage1", "FDR.stage2", "power.stage2", "AUC.stage2",
                    "pAUC.stage2", "FDR.Anno", "power.Anno", "AUC.Anno", "pAUC.Anno")

##### 5. Comparison between LSMM and GPA #####

M      <- 100000   # No. of SNPs
L      <- 10       # No. of fixed effects
K      <- 500      # No. of random effects
alpha  <- 0.2      # parameter in the Beta distribution
Z.perc <- 0.1      # frequency of fixed effects Pr(Z=1)
A.perc <- 0.1      # frequency of random effects Pr(A=1)
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # coefficients of fixed effects
omega  <- 0.2      # proportion of relevant random effects
sigma2 <- 1        # parameter in the spike-slab prior
corr   <- 0.2      # correlation coeffecient among the first 10 random effects
rep    <- 50       # repeat times

result <- matrix(0, rep, 16)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  data <- generate_data_corr(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, corr)

  # LSMM
  fit <- LSMM(data$Pvalue, data$Z, data$A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
  result[i, 1:4]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))

  # LSMM without fixed effects
  fit1 <- LSMM(data$Pvalue, NULL, data$A)
  assoc1.SNP <- assoc.SNP(fit1, FDRset = 0.1, fdrControl = "global")
  result[i, 5:8]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))

  # GPA
  fit.GPA <- GPA(data$Pvalue, data$A)
  assoc.GPA <- assoc(fit.GPA, FDR = 0.1, fdrControl="global")
  out[i, 9:12] <- as.numeric(performance(data$gamma, assoc.GPA, fdr(fit.GPA)))
}

result5 <- as.data.frame(result)
names(result5) <- c("FDR.LSMM", "power.LSMM", "AUC.LSMM", "pAUC.LSMM", "FDR.LSRM", "power.LSRM",
                    "AUC.LSRM", "pAUC.LSRM", "FDR.GPA", "power.GPA", "AUC.GPA", "pAUC.GPA")

##### 6. Comparison between LSMM and cmfdr #####

M      <- 5000     # No. of SNPs
L      <- 5        # No. of fixed effects
alpha  <- 0.2      # parameter in the Beta distribution
Z.perc <- 0.1      # frequency of fixed effects Pr(Z=1)
A.perc <- 0.1      # frequency of random effects Pr(A=1)
beta0  <- -2       # intercept of the logistic model
set.seed(1)
b      <- rnorm(L) # coefficients of fixed effects
omega  <- 0.2      # proportion of relevant random effects
sigma2 <- 1        # parameter in the spike-slab prior
Ks     <- 5        # No. of random effects
rep    <- 10       # repeat times

# 6.1 cmfdr
library(locfdr)
library(pgnorm)
source("run_cmlocfdr.R")

bases_X = NULL;
K = 2
nIter = 2000;
thin = 10;
burnIn = 500;
SSA = 1
SSG = 4
MA = 8
MG = 5
theoNULL = FALSE
mu = .68
inits = NULL

result <- matrix(0, rep, 4)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  data <- generate_data(M, L, Ks, alpha, Z.perc, A.perc, beta0, b, omega, sigma2)
  X <- cbind(rep(1, M), data$Z, data$A)
  colnames(X) <- c("Intcpt", paste("Z", 1:L, sep = ""), paste("A", 1:Ks, sep = ""))

  results=run_cmlocfdr(Pvalue = 1, P = data$Pvalue, X = X,bases_X = bases_X, K = K,
                       nIter = nIter, thin = thin, burnIn = burnIn, SSA = SSA, SSG = SSG,
                       MA = MA, MG = MG, theoNULL = theoNULL, mu = mu, inits = inits)

  ALPHA_array=results[[1]]
  BETA_array=results[[2]]
  GAMMA_array=results[[3]]
  SIGMA_SQ_array=results[[4]]

  first=1
  Alpha_vec=ALPHA_array[[first]]
  Beta_vec=BETA_array[[first]]
  Gamma_vec=GAMMA_array[[first]]
  Alpha=ALPHA_array[[first]]
  Beta=BETA_array[[first]]
  Gamma=GAMMA_array[[first]]
  if(theoNULL==FALSE){Sigma_sq=SIGMA_SQ_array[[first]]}
  for(k in (first+1):length(ALPHA_array)){
    Alpha_vec=cbind(Alpha_vec,ALPHA_array[[k]])
    Beta_vec=cbind(Beta_vec,BETA_array[[k]])
    Gamma_vec=cbind(Gamma_vec,GAMMA_array[[k]])
    Alpha=Alpha+ALPHA_array[[k]]
    Beta=Beta+BETA_array[[k]]
    Gamma=Gamma+GAMMA_array[[k]]
    if(theoNULL==FALSE){Sigma_sq=Sigma_sq+SIGMA_SQ_array[[k]]}
  }
  Alpha=cbind(Alpha/length(first:length(ALPHA_array)))
  Beta=Beta/length(first:length(ALPHA_array))
  Gamma=cbind(Gamma/length(first:length(ALPHA_array)))
  if(theoNULL==FALSE){Sigma_sq=Sigma_sq/length(first:length(ALPHA_array))}
  if(theoNULL==TRUE){Sigma_sq=1}

  Z <- qnorm(data$Pvalue/2)
  f0<-2*dnorm(abs(Z),mean=0,sd=Sigma_sq^.5)
  f1<-f0
  for (l in 1:length(Z)){
    X_l=rbind(X[l,])
    f1[l]<-2*dgamma(abs(Z[l])-mu,shape=exp(X_l%*%Alpha),rate=Beta)
  }
  p1=exp(X%*%Gamma);p0=1;
  pi0=p0/(p0+p1)
  pi1=1-pi0
  cmfdr=pi0*f0/(pi0*f0+pi1*f1)

  est <- rep(0, length(cmfdr))
  cmFDR <- post2FDR(1-cmfdr)
  est[which(cmFDR <= 0.1)] <- 1

  result[i, 1:4] <- as.numeric(performance(data$gamma, est, cmfdr))

}

result6.1 <- as.data.frame(result)
names(result6.1) <- c("FDR", "power", "AUC", "pAUC")

# 6.2 LSMM
result <- matrix(0, rep, 4)

for (i in 1:rep){
  cat(i, "out of", rep, "\n")

  data <- generate_data(M, L, Ks, alpha, Z.perc, A.perc, beta0, b, omega, sigma2)

  fit <- LSMM(data$Pvalue, data$Z, data$A)
  assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")

  result[i, 1:4]  <- as.numeric(performance(data$gamma, assoc.SNP$gamma, 1-fit$pi1))
}

result6.2 <- as.data.frame(result)
names(result6.2) <- c("FDR", "power", "AUC", "pAUC")
