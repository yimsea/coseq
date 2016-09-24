#' Summarize results from clustering using a Normal mixture model
#' 
#' A function to summarize the clustering results obtained from a Normal
#' mixture model estimated using \code{NormMixClus}.
#' 
#' The summary function for an object of class \code{"NormMixClus"}
#' provides the number of clusters selected for the ICL
#' model selection approach.
#' 
#' @param object An object of class \code{"NormMixClus"} 
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables
#' @param digits Integer indicating the number of decimal places to be used
#' for mixture model parameters
#' @param ... Additional arguments
#' @author Andrea Rau
#' @seealso \code{\link{NormMixClus}}, \code{\link{NormMixClus_K}}
#' @references 
#' 
#' Rau, A. and Maugis-Rabusseau, C. (2016) Transformation and model choice for
#' co-expression analayis of RNA-seq data. bioRxiv, doi: http://dx.doi.org/10.1101/065607.
#' @keywords methods
#' @example /inst/examples/NormMixClus.R
#' @export

summary.NormMixClus <-
  function (object, y_profiles, digits=3, ...) 
  {
    x <- object
    if (class(x) != "NormMixClus") {
      stop(paste(sQuote("object"), sep = ""), " must be of class ", 
           paste(dQuote("NormMixClus"), sep = ""), sep = "")
    }
    cat("*************************************************\n")
    clustNum <-  paste(x$nbCluster.all, collapse=",")
    clustErr <- paste(x$nbCluster.error, collapse=",")
    
    
    # clustErr <- paste(rep(paste(x$nbCluster.error, collapse=","),10), collapse=",")
    # strwrap(clustErr, width=5, exdent=2, simplify=FALSE)
                        
                        
                        
    clustErr <- ifelse(clustErr == "", "---", clustErr)
    cat("Clusters fit: ", clustNum, "\n", sep = "")
    cat("Clusters with errors: ", clustErr, "\n", sep= "")
    cat("Selected number of clusters via ICL: ", x$ICL.results$nbCluster, "\n", sep = "")
    cat("ICL of selected model: ", x$ICL.results$ICL, "\n", sep = "")
    cat("*************************************************\n")
    
    x <- object$ICL.results
    
    probaPost <- x$probaPost
    labels <- apply(probaPost, 1, which.max)
    

    
    map <- apply(probaPost, 1, max)
    length(which(map > 0.9))/length(map)
    
    tab <- table(labels)
    names(tab) <- paste("Cluster", names(tab))
    param <- NormMixParam(x=x, y_profiles=y_profiles, digits=digits) 
    mu <- param$mu
    pi <- param$pi
    rownames(mu) <- names(pi) <- names(tab)
    colnames(mu) <- colnames(y_profiles)
    g <- x$nbCluster
    
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
    
    
    
    cat("Mu:\n"); print(round(mu,digits=digits)); cat("\n")
    cat("Pi:\n"); print(round(pi,digits=digits)); cat("\n")
    
    
  }




