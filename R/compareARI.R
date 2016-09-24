#' Pairwise comparisons of ARI values among a set of clustering partitions
#' 
#' Provides the adjusted rand index (ARI) between pairs of clustering paritions.
#'
#' @param x Object of class \code{coseq}, \code{NormMixClus}, or 
#' \code{PoisMixClusWrapper} (object defined in the HTSCluster package that
#' is indirectly called by \code{coseq} for Poisson mixture models),
#' or alternatively a \emph{n} x \emph{M} \code{data.frame} or \code{matrix} 
#' containing the clustering partitions for \emph{M} different models 
#' @param K If \code{NULL}, pairwise ARI values will be calculated among every model in
#' object \code{x}. Otherwise, \code{K} provides a vector of cluster numbers identifying
#' a subset of models in \code{x}.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel 
#' execution using BiocParallel (see next argument \code{BPPARAM}). 
#' Note that parallelization is unlikely to be helpful unless the number of 
#' observations \emph{n} in the clustering partitions or the number of 
#' models \emph{M} are very large.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} 
#' when \code{parallel=TRUE}. If not specified, the parameters last registered 
#' with \code{register} will be used. 
#' @param plot If \code{TRUE}, provide a heatmap using corrplot to visualize 
#' the calculated pairwise ARI values.
#' @param ... Additional optional parameters for corrplot
#'
#' @return Matrix of adjusted rand index values calculated between each pair 
#' of models.
#' 
#' @author Andrea Rau
#' 
#' @export
#' @importFrom HTSCluster highDimensionARI
#' @importFrom corrplot corrplot
#' 
#' @example inst/examples/coseq-package.R
#'
compareARI <- function(x, K=NULL, parallel=FALSE, BPPARAM=bpparam(), plot=TRUE, ...) {
  
  arg.user <- list(...)
  if(is.null(arg.user$digits)) arg.user$digits<-2;
  
  ## For class coseq
  if(class(x) == "coseq") {
    full_labels <- do.call("cbind", lapply(x$results$all.results, 
                                           function(y) apply(y$probaPost,1,which.max)))
  }
  
  ## For class NormMixClus and PoisMixClusWrapper
  if(class(x) == "NormMixClus" | class(x) == "PoisMixClusWrapper") {
    full_labels <- do.call("cbind", lapply(x$all.results, 
                                           function(y) apply(y$probaPost,1,which.max)))
  }
  
  if(!is.null(K)) {
    if(length(K) <= 1) {
      stop("K must be a vector of length at least 2.")
    }
    index <- which(substr(colnames(full_labels), 3, 10) %in% K)
    if(length(index) == 0) stop("None of the indicated models are included in argument x");
    full_labels <- full_labels[,index]
  }
  
  ## For class data.frame or matrix
  if(class(x) == "data.frame" | class(x) == "matrix") {
    if(!is.null(K)) message("K argument only used for objects of class coseq, NormMixClus, and PoisMixClusWrapper.")
    full_labels <- x
    if(length(colnames(full_labels)) == 0) {
      colnames(full_labels) <- paste("Model", 1:ncol(full_labels));
    }
  }
  

  ARI <- matrix(0, nrow=ncol(full_labels), ncol=ncol(full_labels))
  rownames(ARI) <- colnames(ARI) <- colnames(full_labels)

  index <- ARI
  index[upper.tri(index)] <- seq(1, (ncol(ARI) * nrow(ARI) - ncol(ARI))/2)
  
  if(!parallel) {
    tmp <- lapply(1:max(index, na.rm=TRUE), function(ii) {
      index2 <- which(index == ii, arr.ind=TRUE)
      ARItmp <- highDimensionARI(full_labels[,index2[1]], full_labels[,index2[2]])
      return(list(ARItmp=ARItmp, ii=ii))
    })  
  } else if(parallel) {
    tmp <- bplapply(1:max(index, na.rm=TRUE), function(ii) {
      index2 <- which(index == ii, arr.ind=TRUE)
      ARItmp <- highDimensionARI(full_labels[,index2[1]], full_labels[,index2[2]])
      return(list(ARItmp=ARItmp, ii=ii))
    }, BPPARAM=BPPARAM) 
  }
  
  for(i in 1:length(tmp)) {
    new_index <- which(index == tmp[[i]]$ii, arr.ind=TRUE)
    ARI[new_index] <- tmp[[i]]$ARItmp
  }

  diag(ARI) <- 1
  
  if(plot == TRUE) {
    corrplot(ARI, is.corr=FALSE, method="color", type="upper", p.mat=ARI, insig="p-value", 
             sig.level=-1, tl.pos="d", addgrid.col="white",
             tl.col="black", ...)
  }
  
  ARI <- round(ARI, digits = arg.user$digits)
  ARI[lower.tri(ARI)] <- ""
  ARI <- data.frame(ARI, check.names=FALSE)
  return(ARI)
}