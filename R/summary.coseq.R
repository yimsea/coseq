#' Summarize results from clustering using a Poisson or Gaussian mixture model
#' 
#' A function to summarize the clustering results obtained from a Poisson or 
#' Gaussian mixture model estimated using \code{coseq}.
#' 
#' 
#' @param object An object of class \code{"coseq"}
#' @param ... Additional arguments
#' @author Andrea Rau
#' @seealso \code{\link{coseq}}
#' @references 
#' Rau, A. and Maugis-Rabusseau, C. (2016) Transformation and model choice for
#' co-expression analayis of RNA-seq data. bioRxiv, doi: http://dx.doi.org/10.1101/065607.
#' 
#' @keywords methods
#' @example /inst/examples/coseq-package.R
#' @export
summary.coseq <-
  function (object,  ...) 
  {
    x <- object
    if (class(x) != "coseq") {
      stop(paste(sQuote("object"), sep = ""), " must be of class ", 
           paste(dQuote("coseq"), sep = ""), sep = "")
    }
    cat("*************************************************\n")
    cat("Model: ", x$model, "\n", sep = "")
    cat("Transformation: ", x$transformation, "\n", sep = "")

    if(class(x$results) == "NormMixClus") summary(object=x$results, y_profiles=x$tcounts, ...)
    if(class(x$results) == "PoisMixClus") summary(x$results, ...)
  }




