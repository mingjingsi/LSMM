\name{LSMM}
\alias{LSMM}
\title{
  Fit LSMM
}
\description{
Fit LSMM.
}
\usage{
LSMM(Pvalue, Z = NULL, A = NULL, alpha = 0.1, pi1_ = 0.1, maxiter = 10000, tol = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pvalue}{
     p-value from GWAS data, the length is the number of SNPs.
}
  \item{Z}{
     Design matrix of fixed effects without intercept, where row and column correspond to SNP and annotation, respectively. Default is \code{NULL}.
}
  \item{A}{
     Design matrix of random effects, where row and column correspond to SNP and annotation, respectively. Default is \code{NULL}.
}
  \item{alpha}{
     Initial value of parameter \code{alpha} in the Beta distribution which is between 0 and 1. Default is 0.1.
}
  \item{pi1_}{
     Initial value of parameter pi1 which is the proportion that gamma equals 1. Default is 0.1.
}
  \item{maxiter}{
     Maximum number of EM iteration. Default is 10000.
}
  \item{tol}{
     Threshold to stop the EM iteration. Default is 1e-06.
}
}
\details{
  \code{LSMM} fits the LSMM. It requires to provide GWAS p-value to \code{Pvalue}, while users can also provide annotation data to \code{Z} and \code{A} which correspond to fixed-effects and random-effects respectively. It is assumed that the length of \code{Pvalue} equals to the number of rows of matrix provided to \code{Z} and \code{A}.

  \code{assoc.SNP} makes statistical inference of risk SNPs. It requires LSMM model fit. The results can be obtained by setting the type of FDR and the corresponding FDR threshold.

  \code{relev.Anno} makes statistical inference of relevant tissue-specific functional annotations. It requires LSMM model fit. The results can be obtained by setting the type of FDR and the corresponding FDR threshold.
}
\value{
List of parameter estimation, the value of lower bound of Log-likelihood and iteration times of each stage.
  \item{alpha}{
     Estimation of parameter \code{alpha} in the Beta distribution.
}
  \item{pi1_}{
     Estimation of parameter pi1 in Two Groups Model, which measures the proportion that gamma equals 1.
}
  \item{pi1.stage1}{
     Posterior probability in the stage1 (i.e., under the Two Groups Model) that each SNP is associated with this phenotype.
}
  \item{pi1.stage2}{
     Posterior probability in the stage2 (i.e., under the LFM) that each SNP is associated with this phenotype.
}
  \item{pi1}{
     Posterior probability in the main stage that each SNP is associated with this phenotype.
}
  \item{b}{
     Estimation of parameter \code{b} which is the coefficients of fixed effects.
}
  \item{sigma2}{
     Estimation of parameter \code{sigma2} which is the parameter in spike-slab prior.
}
  \item{omega}{
     Estimation of parameter \code{omega} which measures the proportion of relevant annotations
corresponding to this phenotype.
}
  \item{omegak}{
     Posterior probability that each annotation is relevant to this phenotype.
}
  \item{beta}{
     Estimation of effect size of each annotation.
}
  \item{Lq}{
     The value of lower bound of Log-likelihood in the main stage.
}
  \item{iter_times.stage1}{
     Iteration times in the stage1.
}
  \item{iter_times.stage2}{
     Iteration times in the stage2.
}
  \item{iter_times.stage3}{
     Iteration times in the stage3.
}
  \item{iter_times.stage4}{
     Iteration times in the stage4.
}
  \item{iter_times}{
     Iteration times in the main stage under the Two Groups Model.
}


}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming and Can Yang
}
\seealso{
 \code{\link{assoc.SNP}}, \code{\link{relev.Anno}}
}
\examples{
data(ExampleData)

Pvalue <- ExampleData$Pvalue
Z <- ExampleData$Z
A <- ExampleData$A

fit <- LSMM(Pvalue, Z, A)

}
