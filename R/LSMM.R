# main LSMM function

# Inputs:
# Pvalue : a M by 1 vector of p-values of the phenotype, M is the number of SNPs
# Z      : a M by L matrix, fixed effects (not include intercept), L is the number of fixed effects, for example the genic annotation categories, set to NULL by default
# A      : a M by K matrix, random effects, K is the number of random effects, for example the tissue-specific functional annotations, set to NULL by default
# alpha  : parameter in Beta(alpha, 1), set initial value to 0.1 by default
# pi1_   : the proportion of risk SNPs, set initial value to 0.1 by default
# maxiter: maximize iterations, set to 1e4 by default
# tol    : tolerance, set to 1e-6 by default

# Outputs:
# parameter estimations: alpha, pi1_(no Z and no A), b, sigma2, omega, beta
# posterior of latent variables: pi1, omegak
# Lower bound evaluation: Lq
# iteration times: iter_times

LSMM <- function(Pvalue, Z=NULL, A=NULL, alpha=0.1, pi1_=0.1, maxiter=1e4, tol=1e-6){

  # check whether p-values are in [0, 1]
  if (any(Pvalue < 0 | Pvalue > 1)) {
    stop("Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one." )
  }

  M <- length(Pvalue)

  # check whether fixed effects are inputted, convert Z to matrix and add a colume of 1 as the intercept
  if (is.null(Z)) {
    L <- 0
  }
  else {
    if (!is.matrix(Z)) {
      Z <- as.matrix(Z)
    }
    L <- ncol(Z)
  }
  Z <- cbind(rep(1, M), Z)

  # check whether random effects (annotations) are inputted and convert A to matrix
  if (is.null(A)) {
    K <- 0
  }
  else {
    if (!is.matrix(A)) {
      A <- as.matrix(A)
    }
    K <- ncol(A)
  }

  # report setting
  message("Info: Number of SNPs: ", M)
  message("Info: Number of fixed effects: ",  L)
  message("Info: Number of random effects: ", K)

  # set zero p-values to small values to avoid log(0)
  if (any(Pvalue < 1e-30)) {
    message("Info: Some SNPs have p-values close to zero." )
    message("Info: Number of SNPs with p-values close to zero: ", length(which(Pvalue < 1e-30)))
    message("Info: p-values for these SNPs are set to ", 1e-30 )

    Pvalue[Pvalue < 1e-30] <- 1e-30
  }

  # fit LSMM
  if (is.null(A)) {
    if (L == 0) {
      # If no fixed effects and random effects, use Two-Groups Model
      fit <- TwoGroup_Rcpp(Pvalue, alpha, pi1_, maxiter, tol)
    }
    else
      # If no random effects, use Latent Fixed-effect Model
      fit <- LFM_Rcpp(Pvalue, Z, alpha, pi1_, maxiter, tol)
  }
  else
    fit <- LSMM_Rcpp(Pvalue, Z, A, alpha, pi1_, maxiter, tol)

  return(fit)

}
