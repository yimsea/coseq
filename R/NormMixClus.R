#' Normal mixture model estimation and selection for a series of cluster numbers
#'
#' Perform co-expression and co-abudance analysis of high-throughput 
#' sequencing data, with or without data transformation, using a Normal 
#' mixture models. The output of \code{NormMixClus} is an S3 object of 
#' class \code{NormMixClus}.
#'
#' @param y_profiles  (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value or a sequence of values).
#' @param subset Optional vector providing the indices of a subset of
#' genes that should be used for the co-expression analysis (i.e., row indices
#' of the data matrix \code{y}.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel 
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running 
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects 
#' from the current R environment before calling the function, as it is possible that R's 
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when 
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register} 
#' will be used.
#' @param ... Additional optional parameters to be passed to \code{\link{NormMixClus_K}}.
#'
#' @return 
#' An S3 object of class \code{NormMixClus} containing the following:
#' \item{nbCluster.all }{Vector giving the number of clusters for each of the fitted models}
#' \item{loglike.all }{Log likelihoods calculated for each of the fitted models}
#' \item{ICL.all}{ICL values calculated for each of the fitted models}
#' \item{ICL.results }{Object of class
#' \code{NormMixClus} giving the results from the model chosen via the ICL
#' criterion} 
#' \item{all.results }{List of objects of class \code{NormMixClus} giving the results for all models}
#' 
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#' 
#' @example inst/examples/NormMixClus.R
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @export

NormMixClus <- function(y_profiles, K, subset=NULL, parallel=TRUE, BPPARAM=bpparam(), ...){
  
  subset.index <- subset
  
  ## Parse ellipsis function
  arg.user <- list(...)
  if(is.null(arg.user$alg.type)) arg.user$alg.type<-"EM";
  if(is.null(arg.user$init.runs)) arg.user$init.runs<-50;
  if(is.null(arg.user$init.type)) arg.user$init.type<-"small-em";
  if(is.null(arg.user$init.iter)) arg.user$init.iter<-20;
  if(is.null(arg.user$iter)) arg.user$iter<-1000;
  if(is.null(arg.user$cutoff)) arg.user$cutoff<-0.001;
  if(is.null(arg.user$GaussianModel)) arg.user$GaussianModel<-"Gaussian_pk_Lk_Ck";
  if(is.null(arg.user$verbose)) arg.user$verbose<-TRUE;
  if(is.null(arg.user$digits)) arg.user$digits <- 3;
  
  y_profiles <- as.data.frame(y_profiles)
  
  ## In case only a subset of data are to be used for analysis
  if(!is.null(subset.index)) {
    y_profiles <- y_profiles[subset.index,]
  }
                              
  all.results <- vector("list", length = length(K))
  names(all.results) <- paste("K=", K, sep = "")
  if(arg.user$verbose == TRUE) {
    cat("Running K =", min(K), "...\n")
  }
  all.results[[1]] <- suppressWarnings(NormMixClus_K(y_profiles, K=min(K), alg.type=arg.user$alg.type, 
                                                     init.runs=arg.user$init.runs,
                                                     init.type=arg.user$init.type, 
                                                     init.iter=arg.user$init.iter,
                                                     iter=arg.user$iter, 
                                                     cutoff=arg.user$cutoff,
                                                     GaussianModel=arg.user$GaussianModel, 
                                                     digits=arg.user$digits))
  index <- 2
  remainingK <- K[-which(K == min(K))]
  if (length(remainingK) > 0) {
    ## In the case where parallelization is NOT used
    if(!parallel) {
      for (k in remainingK) {
        if(arg.user$verbose == TRUE) {
          cat("Running K =", k, "...\n")
        }
        all.results[[index]] <- suppressWarnings(NormMixClus_K(y_profiles=y_profiles, K=k, 
                                                               alg.type=arg.user$alg.type, 
                                                               init.runs=arg.user$init.runs,
                                                               init.type=arg.user$init.type, 
                                                               init.iter=arg.user$init.iter,
                                                               iter=arg.user$iter, 
                                                               cutoff=arg.user$cutoff,
                                                               GaussianModel=arg.user$GaussianModel, 
                                                               digits=arg.user$digits))
        index <- index + 1
      }
      ## In the case where parallelization IS used
    } else if(parallel) {
      tmp <- bplapply(remainingK, function(ii) {
        if(arg.user$verbose == TRUE) {
          cat("Running K =", ii, "...\n")
        }
        res <- suppressWarnings(NormMixClus_K(y_profiles=y_profiles, K=as.numeric(ii),
                             alg.type=arg.user$alg.type, init.runs=arg.user$init.runs,
                             init.type=arg.user$init.type, init.iter=arg.user$init.iter,
                             iter=arg.user$iter, cutoff=arg.user$cutoff,
                             GaussianModel=arg.user$GaussianModel, digits=arg.user$digits))

        return(res)}, BPPARAM=BPPARAM)
      Kmods <- paste("K=", unlist(lapply(tmp, function(x) x$nbCluster)), sep="")
      # cat(Kmods, "\n")
      # cat(length(tmp), "\n")
      # cat(length(match(Kmods, names(all.results))), "\n")
      # print(unlist(lapply(tmp, function(xx) length(xx$K) != 0)))
      all.results[match(Kmods, names(all.results))] <- tmp[which(unlist(lapply(tmp, function(xx) 
        length(xx$nbCluster) != 0)) == TRUE)]
    }
  }
  
  nbClust.all<-unlist(lapply(all.results, function(x) x$nbCluster))
  logLike.all <- unlist(lapply(all.results, function(x) x$log.like))
  ICL.all <- unlist(lapply(all.results, function(x) x$ICL))
  ICL.choose <- names(ICL.all)[which.min(ICL.all)]
  select.results <- all.results[[ICL.choose]]
  
  all.results <- lapply(all.results, function(xx) {
    if(is.null(xx) == TRUE) {
        xx <- list(probaPost = matrix(NA,nrow=0,ncol=0), 
               log.like = vector("numeric", 0), ICL = vector("numeric", 0), 
               K = vector("numeric", 0))
    class(xx) <- "NormMixClus_K"
    }
    return(xx)
  })
  
  RES <- list(nbCluster.all=nbClust.all,logLike.all = logLike.all, ICL.all = ICL.all, 
              all.results = all.results, ICL.results = select.results,
              nbCluster.error=K[!K %in% nbClust.all])
  class(RES) <- "NormMixClus"
  return(RES)
}