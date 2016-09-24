#' Normal mixture model estimation
#'
#' Perform co-expression and co-abudance analysis of high-throughput 
#' sequencing data, with or without data transformation, using a Normal 
#' mixture models for single number of clusters \emph{K}. 
#' The output of \code{NormMixClus_K} is an S3 object of 
#' class \code{NormMixClus_K}.
#'
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value).
#' @param init.type Type of initialization strategy to be used:
#' \dQuote{\code{small-em}} for the Small-EM strategy, \dQuote{\code{random}}, \dQuote{\code{CEM}},
#' or \dQuote{\code{SEMMax}}
#' @param GaussianModel One of the 28 forms of Gaussian models defined in Rmixmod,
#' by default equal to the \code{"Gaussian_pk_Lk_Ck"} (i.e., a general family model with free
#' proportions, free volume, free shape, and free orientation)
#' @param init.runs Number of runs to be used for the Small-EM strategy, with a default value of 50
#' @param init.iter Number of iterations to be used within each run for the
#' Small-EM strategry, with a default value of 20
#' @param alg.type Algorithm to be used for parameter estimation:
#' \dQuote{\code{EM}}, \dQuote{\code{CEM}}, \dQuote{\code{SEM}}
#' @param cutoff Cutoff to declare algorithm convergence
#' @param iter Maximum number of iterations to be run for the chosen algorithm
#' @param verbose If \code{TRUE}, verbose output is created
#' @param digits Integer indicating the number of decimal places to be used for the
#' \code{probaPost} output
#'
#' @return
#' An S3 object of class \code{NormMixClus_K} containing the following:
#' \item{probaPost }{Matrix containing the conditional
#' probabilities of belonging to each cluster for all observations}
#' \item{log.like }{Value of log likelihood} 
#' \item{ICL }{Value of ICL criterion} 
#' \item{nbCluster }{Number of clusters (equivalent to \code{K})}
#' \item{GaussianModel }{Gaussian model form fit in the mixture model}
#' 
#' @author Cathy Maugis-Rabusseau, Andrea Rau
#' 
#' @importFrom Rmixmod mixmodGaussianModel
#' @importFrom Rmixmod mixmodStrategy
#' @importFrom Rmixmod mixmodCluster
#' @example inst/examples/NormMixClus.R
#' 
#' @export

NormMixClus_K <- function(y_profiles, K, alg.type="EM", init.runs=50, 
                          init.type="small-em", GaussianModel="Gaussian_pk_Lk_Ck",
                          init.iter=20, iter=1000, cutoff=0.001, verbose=TRUE, digits=3) {
     
  if(!is.data.frame(y_profiles)) y_profiles <- as.data.frame(y_profiles)
  
  models <- mixmodGaussianModel(listModels=c(GaussianModel))
  # strategy   
  if(init.type == "small-em") init.type <- "smallEM";
  strats<-mixmodStrategy(algo=alg.type, nbTry = 1,
                         initMethod = init.type, nbTryInInit = init.runs,
                         nbIterationInInit = init.iter, nbIterationInAlgo = iter,
                         epsilonInInit = cutoff, epsilonInAlgo = cutoff,
                         seed = NULL)
  xem <- mixmodCluster(y_profiles, nbCluster=K, models=models, strategy=strats, 
                                      criterion="ICL")
  #param<-list(pi=xem["bestResult"]@parameters@proportions,
  #            mu=xem["bestResult"]@parameters@mean,
  #            variance=xem["bestResult"]@parameters@variance)
  pp <- round(xem["bestResult"]@proba, digits=digits )
  
  if(nrow(pp) != 0 & !is.null(rownames(y_profiles))) rownames(pp) <- rownames(y_profiles)
  
  res <- list(probaPost=pp, log.like=xem["bestResult"]@likelihood,
              ICL=xem["bestResult"]@criterionValue, nbCluster=xem["bestResult"]@nbCluster,
              GaussianModel=GaussianModel)
  class(res) <- "NormMixClus_K"
  return(res)
}
