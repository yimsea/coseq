#' Summarize results from clustering using a Normal mixture model
#' 
#' A function to summarize the clustering results obtained from a Normal
#' mixture model.
#' 
#' The summary function for an object of class \code{"NormMixClus_K"} provides the
#' following summary of results:
#' 
#' 1) Number of clusters and model selection criterion used, if applicable.
#' 
#' 2) Number of observations across all clusters with a maximum conditional
#' probability greater than 90% (and corresponding percentage of total
#' observations) for the selected model.
#' 
#' 3) Number of observations per cluster with a maximum conditional probability
#' greater than 90% (and corresponding percentage of total observations per
#' cluster) for the selected model.
#' 
#' 4) \eqn{\ensuremath\boldsymbol{\mu}}{\mu} values for the selected
#' model.
#' 
#' 5) \eqn{\ensuremath\boldsymbol{\pi}}{\pi} values for the selected model.
#' 
#' @param object An object of class \code{"NormMixClus_K"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables
#' @param digits Integer indicating the number of decimal places to be used
#' for mixture model parameters
#' @param ... Additional arguments
#' 
#' @author Andrea Rau
#' @seealso \code{\link{NormMixClus}}, \code{\link{NormMixClus_K}}
#' @keywords methods
#' @example /inst/examples/NormMixClus.R
#' @export
`summary.NormMixClus_K` <-
  function (object, y_profiles, digits=3, ...) 
  {
    x <- object
    if (class(x) != "NormMixClus_K") {
      stop(paste(sQuote("object"), sep = ""), " must be of class ", 
           paste(dQuote("NormMixClus_K"), sep = ""), sep = "")
    }
    
    probaPost <- x$probaPost
    labels <- apply(probaPost, 1, which.max)
    
    param <- NormMixParam(x, y_profiles, digits=digits) 
    mu <- param$mu
    pi <- param$pi
    g <- x$nbCluster
    
    map <- apply(probaPost, 1, max)
    length(which(map > 0.9))/length(map)
    
    cat("*************************************************\n")
    cat("Number of clusters = ", g, "\n", sep = "")
    cat("ICL = ", x$ICL, "\n", sep = "")
    cat("*************************************************\n")
    tab <- table(labels)
    names(tab) <- paste("Cluster", names(tab))
    cat("Cluster sizes:\n"); print(tab); cat("\n")
    cat("Number of observations with MAP > 0.90 (% of total):\n")
    cat(length(which(map > 0.9)), " (", round(length(which(map > 0.9))/length(map)*100,2),
        "%)\n\n", sep = "")
    cat("Number of observations with MAP > 0.90 per cluster (% of total per cluster):\n"); 
    
    tab2 <- matrix(NA, nrow = 2, ncol = g)
    colnames(tab2) <- paste("Cluster", 1:g); rownames(tab2) <- rep("", 2)
    for(i in 1:g) {
      if(sum(labels == i) > 1) {
        map.clust <- apply(matrix(probaPost[labels == i,], ncol=g), 1, max)
        tab2[1,i] <- length(which(map.clust > 0.9))
        tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
                           "%)", sep = "")
      }
      if(sum(labels == i) == 1) {
        map.clust <- max(probaPost[labels == i,])
        tab2[1,i] <- length(which(map.clust > 0.9))
        tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
                           "%)", sep = "")
      }
      if(sum(labels == i) == 0) {
        tab2[1,i] <- "---"
        tab2[2,i] <- "---"
      }
    }
    print(tab2, quote = FALSE); cat("\n")
    
    rownames(mu) <- names(pi) <- names(tab)
    colnames(mu) <- colnames(y_profiles)
    
    
    cat("Mu:\n"); print(round(mu,digits=digits)); cat("\n")
    cat("Pi:\n"); print(round(pi,digits=digits)); cat("\n")
  }
