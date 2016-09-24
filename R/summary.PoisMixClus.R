#' Summarize results from clustering using a Poisson mixture model
#' 
#' A function to summarize the clustering results obtained from a Poisson
#' mixture model estimated using \code{coseq} (which indirectly calls
#' \code{HTSCluster}.
#' 
#' The summary function for an object of class \code{"PoisMixClus"}
#' provides the number of clusters selected according to the user-defined
#' model selection criterion.
#' 
#' @param object An object of class \code{"PoisMixClus"} 
#' @param digits Integer indicating the number of decimal places to be used
#' for mixture model parameters
#' @param ... Additional arguments
#' @author Andrea Rau
#' @seealso \code{\link{coseq}}
#' @references Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux,
#' G. (2015) Co-expression analysis of high-throughput transcriptome sequencing
#' data with Poisson mixture models. Bioinformatics, doi:
#' 10.1093/bioinformatics/btu845.
#' 
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011).
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' 
#' @keywords methods
#' @example inst/examples/PoisMixClus_ex.R
#' @export

summary.PoisMixClus <-
  function (object, digits=3, ...) 
  {
    x <- object
    if (class(x) != "PoisMixClus") {
      stop(paste(sQuote("object"), sep = ""), " must be of class ", 
           paste(dQuote("PoisMixClus"), sep = ""), sep = "")
    }
    cat("*************************************************\n")
    clustNum <-  paste(x$nbCluster.all, collapse=",")
    
    cat("Clusters fit: ", clustNum, "\n", sep = "")
    cat("Selected number of clusters: ", 
        ncol(x$selected.results$probaPost), "\n", sep = "")
    cat("Model selection criterion: ", 
        x$selected.results$model.selection, "\n", sep = "")
    cat("*************************************************\n")
    
    x <- object$selected.results
    
    probaPost <- x$probaPost
    labels <- apply(probaPost, 1, which.max)
    
    map <- apply(probaPost, 1, max)
    length(which(map > 0.9))/length(map)
    
    tab <- table(labels)
    names(tab) <- paste("Cluster", names(tab))
    lambda <- x$lambda
    pi <- x$pi
    g <- ncol(lambda)
    
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
    
    cat("lamba:\n"); print(round(lambda,digits=digits)); cat("\n")
    cat("pi:\n"); print(round(pi,digits=digits)); cat("\n")
    
  }




