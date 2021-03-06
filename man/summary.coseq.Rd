% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.coseq.R
\name{summary.coseq}
\alias{summary.coseq}
\title{Summarize results from clustering using a Poisson or Gaussian mixture model}
\usage{
\method{summary}{coseq}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"coseq"}}

\item{...}{Additional arguments}
}
\description{
A function to summarize the clustering results obtained from a Poisson or 
Gaussian mixture model estimated using \code{coseq}.
}
\examples{
## Simulate toy data, n = 300 observations
set.seed(12345)
countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
countmat <- countmat[which(rowSums(countmat) > 0),]
conds <- rep(c("A","B","C","D"), each=2)

## Run the Normal mixture model for K = 2,3,4
run_arcsin <- coseq(y=countmat, K=2:4, iter=5, transformation="arcsin")

## Plot and summarize results
plot(run_arcsin)
summary(run_arcsin)

## Compare ARI values for all models (no plot generated here)
ARI <- compareARI(run_arcsin, plot=FALSE)

## Compare ICL values for models with arcsin and logit transformations
run_logit <- coseq(y=countmat, K=2:4, iter=5, transformation="logit")
compareICL(list(run_arcsin, run_logit))

}
\author{
Andrea Rau
}
\references{
Rau, A. and Maugis-Rabusseau, C. (2016) Transformation and model choice for
co-expression analayis of RNA-seq data. bioRxiv, doi: http://dx.doi.org/10.1101/065607.
}
\seealso{
\code{\link{coseq}}
}
\keyword{methods}

