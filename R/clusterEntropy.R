#' Calculation of per-cluster entropy
#' 
#' Provides the calculation of per-cluster entropy, equivalent to 
#' \deqn{Entropy(k) = \sum_{i \in C_k} \log (\tau_{ik})}
#' where \eqn{\tau_{ik}} is the conditional probability of gene \emph{i} belonging
#' to cluster \emph{k} and \eqn{C_k} corresponds to the set of indices of genes
#' attributed to cluster \emph{k}.
#'
#' @param probaPost Matrix containing the conditional probabilities of belonging 
#' to each cluster for all observations
#'
#' @return Entropy per cluster
#'
#' @author Cathy Maugis-Rabusseau
#' 
#' @examples
#' ## Generate artificial matrix of conditional probabilities for K=5 clusters
#' tmp <- matrix(runif(100*5), nrow=100, ncol=5)
#' probaPost <- tmp / rowSums(tmp)
#' clusterEntropy(probaPost)
#' 
#' @export
clusterEntropy <- function(probaPost) {
  label<-apply(probaPost,1,which.max)
  entrop<-NULL
  for (k in 1:max(label)){
    I <- which(label==k)
    entrop <- c(entrop,sum(log(probaPost[I,k])))
  }
  return(entrop)
}
