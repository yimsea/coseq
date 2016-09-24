#' Calculate the mean and covariance for a Normal mixture model
#' 
#' Calculates the mean and covariance parameters for a normal mixture model
#' of the form pK_Lk_Ck
#'
#' @param x Object of class \code{coseq}, \code{NormMixClus}, or \code{NormMixClus_K}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables, required for \code{x} of class \code{NormMixClus}
#' or \code{NormMixClus_K}
#' @param K The model used for parameter estimation for objects \code{x} of
#' class \code{coseq} or \code{NormMixClus}. When \code{NULL}, the model selected
#' by the ICL criterion is used; otherwise, \code{K} should designate the number
#' of clusters in the desired model
#' @param digits Integer indicating the number of decimal places to be used for output
#' @param plot If \code{true}, produce heatmaps to visualize the estimated per-cluster
#' correlation matrices
#' @param ... Additional optional parameters to pass to \code{corrplot}, if desired
#'
#'
#' @return
#' \item{pi }{ Vector of dimension \emph{K} with the estimated cluster proportions from
#' the Gaussian mixture model, where \emph{K} is the number of clusters}
#' \item{mu }{ Matrix of dimension \emph{K} x \emph{d} containing the estimated mean
#' vector from the Gaussian mixture model, where \emph{d} is the
#' number of samples in the data \code{y_profiles} and \emph{K} is the number of clusters}
#' \item{Sigma }{ Array of dimension \emph{d} x \emph{d} x \emph{K} containing the
#' estimated covariance matrices from the Gaussian mixture model, where \emph{d} is the
#' number of samples in the data \code{y_profiles} and \emph{K} is the number of clusters}
#' \item{rho }{ Array of dimension \emph{d} x \emph{d} x \emph{K} containing the
#' estimated correlation matrices from the Gaussian mixture model, where \emph{d} is the
#' number of samples in the data \code{y_profiles} and \emph{K} is the number of clusters}
#' 
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#' 
#' @example inst/examples/NormMixClus.R
#' @export
#' @importFrom stats cov2cor
#' @importFrom grDevices n2mfrow
#' @importFrom graphics mtext par
#'

NormMixParam <- function(x, y_profiles=NULL, K=NULL, digits=3, plot=FALSE, ...) {
  
  if (class(x) != "coseq" & class(x) != "NormMixClus" & class(x) != "NormMixClus_K") {
    stop(paste(sQuote("x"), sep = ""), " must be of class ", 
         paste(dQuote("coseq"), sep = ""), " or ", paste(dQuote("NormMixClus"), sep = ""),
         " or ", paste(dQuote("NormMixClus_K"), sep = ""), sep = "")
  }
  if(class(x) == "coseq") {
    y_profiles <- x$tcounts
    if(is.character(K) == TRUE) {
      if(K == "ICL") {
        mod <- x$results$ICL.results;
        GaussianModel <- x$results$ICL.results$GaussianModel
      }
    }
    if(is.null(K) == TRUE) {
      mod <- x$results$ICL.results;
      GaussianModel <- x$results$ICL.results$GaussianModel
    }
    if(is.numeric(K) == TRUE) {
      mod <- x$results$all.results[[paste("K=", K,sep="")]]
      GaussianModel <- x$results$all.results[[paste("K=", K,sep="")]]$GaussianModel
    }
    probaPost <- mod$probaPost
  }
  if(class(x) == "NormMixClus") {
    if(is.null(y_profiles) == TRUE) stop(paste(dQuote("y_profiles"), sep = ""), " must
                                           be included for class ", 
                                         paste(dQuote("NormMixClus"), sep = ""))
    if(K == "ICL" | is.null(K)==TRUE) {
      mod <- x$ICL.results;
      GaussianModel <- x$ICL.results$GaussianModel
    }
    if(is.numeric(K) == TRUE) {
      mod <- x$all.results[[paste("K=", K,sep="")]]
      GaussianModel <- x$all.results[[paste("K=", K,sep="")]]$GaussianModel
    }
    probaPost <- mod$probaPost
  }
  if(class(x) == "NormMixClus_K") {
    if(is.null(y_profiles) == TRUE) stop(paste(dQuote("y_profiles"), sep = ""), " must
                                         be included for class ", 
                                         paste(dQuote("NormMixClus_K"), sep = ""))
    mod <- x
    probaPost <- mod$probaPost
    GaussianModel <- x$GaussianModel
  }
  
  if(GaussianModel != "Gaussian_pk_Lk_Ck") 
    stop("Recalculation of Gaussian parameters currently only supported for",
         paste(dQuote("Gaussian_pk_Lk_Ck")))
    
  pi <- apply(probaPost,2,sum) / nrow(y_profiles)
  mu <- matrix(0, nrow=ncol(probaPost), ncol=ncol(y_profiles))
  Sigma <- array(0, dim=c(ncol(y_profiles), ncol(y_profiles), ncol(probaPost)))
  for (k in 1:ncol(probaPost)){
    mu[k,]<-apply(probaPost[,k] * y_profiles,2,sum) / sum(probaPost[,k])
    Sigma[,,k] <- (t(y_profiles) - mu[k,]) %*% ( t(t(y_profiles) - mu[k,]) * probaPost[,k])
    Sigma[,,k] <- Sigma[,,k] / sum(probaPost[,k])
  }
  rho <- lapply(1:ncol(probaPost), function(xx) cov2cor(Sigma[,,xx]))
  rho <- array(unlist(rho), dim = c(dim(rho[[1]]), length(rho)))
  
  dimnames(rho) <- list(colnames(y_profiles), colnames(y_profiles), 
                        paste("Cluster", 1:ncol(probaPost)))
  dimnames(Sigma) <- list(colnames(y_profiles), colnames(y_profiles), 
                        paste("Cluster", 1:ncol(probaPost)))
  colnames(mu) <- colnames(y_profiles)
  rownames(mu) <- paste("Cluster", 1:ncol(probaPost))
  
  if(plot == TRUE) {
    par(mfrow=n2mfrow(ncol(probaPost)))
    for(kk in 1:ncol(probaPost)) {
      corrplot(rho[,,kk],  method="ellipse", type="upper",
               tl.pos="d", ...)
      mtext(paste("K =", kk), side=2, las=1, line=-1)
    }
  }

  param <- list(pi=round(pi, digits), mu=round(mu, digits), Sigma=round(Sigma, digits),
                rho=round(rho, digits))
  return(param)           
}
