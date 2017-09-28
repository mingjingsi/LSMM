\name{LSMM-package}
\alias{LSMM-package}
\docType{package}
\title{
LSMM (Latent Sparse Mixed Model)
}
\description{
This package provides functions for fitting LSMM,
a statistical approach to integrating functional
annotations with genome-wide association studies.
}
\details{
\tabular{ll}{
Package: \tab LSMM\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2017-09-04\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package contains three main function, \code{LSMM}, \code{assoc.SNP} and \code{relev.Anno} . \code{LSMM} fits LSMM, \code{assoc.SNP} makes statistical inference for risk SNPs and \code{relev.Anno} makes statistical inference for relevant tissue-specific functional annotations.
}
\author{
Jingsi Ming

Maintainer: Jingsi Ming <jsming@life.hkbu.edu.hk>
}
\references{

}
\keyword{ package }
\seealso{
\code{\link{LSMM}}, \code{\link{assoc.SNP}}, \code{\link{relev.Anno}}.
}
\examples{
library(LSMM)

data(ExampleData)

Pvalue <- ExampleData$Pvalue
Z <- ExampleData$Z
A <- ExampleData$A

fit <- LSMM(Pvalue, Z, A)
assoc.SNP <- assoc.SNP(fit, FDRset = 0.1, fdrControl = "global")
relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "local")
}
