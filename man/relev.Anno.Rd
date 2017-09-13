\name{relev.Anno}
\alias{relev.Anno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical inference for relevant annotations
}
\description{
Statistical inference for relevant annotations.
}
\usage{
relev.Anno(fit, FDRset = 0.1, fdrControl)
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
\code{relev.Anno} makes statistical inference for relevant tissue-specific functional annotations.
}
\value{
List of binary values indicating relevant annotations for the phenotype, where one indicates that the annotation is relevant to the phenotype and zero otherwise
}
\author{
Jingsi Ming and Can Yang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LSMM}}, \code{\link{assoc.SNP}}
}
\examples{
relev.Anno <- relev.Anno(fit, FDRset = 0.1, fdrControl = "local")

}
