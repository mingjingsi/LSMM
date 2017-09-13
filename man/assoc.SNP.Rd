\name{assoc.SNP}
\alias{assoc.SNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Statistical inference for risk SNPs
}
\description{
  Statistical inference for risk SNPs.
}
\usage{
assoc.SNP(fit, FDRset = 0.1, fdrControl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
     LSMM model fit.
}
  \item{FDRset}{
     FDR level.
}
  \item{fdrControl}{
     Method to control FDR. Possible values are "global" (global FDR control) and "local" (local FDR control).
}
}
\details{
  \code{assoc.SNP} makes statistical inference for risk SNPs.
}
\value{
List of binary values indicating association of SNPs for the phenotype, where one indicates that the SNP is associated with the phenotype and zero otherwise
  \item{gamma.stage1}{
     Estimation of latent status that each SNP is associated with this phenotype in the stage1 (i.e., under the Two Groups Model).
}
  \item{gamma.stage2}{
     Estimation of latent status that each SNP is associated with this phenotype in the stage2 (i.e., under the LFM).
}
  \item{gamma}{
     Estimation of latent status that each SNP is associated with this phenotype in the main stage (i.e., under the LSMM).
}
}
\author{
Jingsi Ming and Can Yang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{LSMM}}, \code{\link{relev.Anno}}
}
\examples{
assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
}
