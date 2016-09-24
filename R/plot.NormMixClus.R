#' Visualize results from clustering using a Normal mixture model
#' 
#' Plot a NormMixClus object.
#'
#' @param x An object of class \code{"NormMixClus"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables to be used for graphing results (optional for
#' \code{logLike}, \code{ICL}, \code{probapost_boxplots}, and \code{probapost_barplots})
#' @param K If desired, the specific model to use for plotting. If \code{NULL},
#' the model chosen by ICL will be plotted
#' @param threshold Threshold used for maximum conditional probability; only observations
#' with maximum conditional probability greater than this threshold are visualized 
#' @param conds Condition labels, if desired
#' @param average_over_conds If \code{TRUE}, average values of \code{y_profiles} within
#' each condition identified by \code{conds} for the \code{profiles} and \code{boxplots}
#' plots
#' @param graphs Graphs to be produced, one (or more) of the following: 
#' \code{"logLike"} (log-likelihood plotted versus number of clusters),
#' \code{"ICL"} (ICL plotted versus number of clusters), 
#' \code{"profiles"} (line plots of profiles in each cluster), \code{"boxplots"} 
#' (boxplots of profiles in each cluster), \code{"probapost_boxplots"} (boxplots of
#' maximum conditional probabilities per cluster), \code{"probapost_barplots"} 
#' (number of observations with a maximum conditional probability greater than 
#' \code{threshold} per cluster), \code{"probapost_histogram"} (histogram of maximum
#' conditional probabilities over all clusters) ...
#' @param order If \code{TRUE}, order clusters in \code{probapost_boxplot} by median and
#' \code{probapost_barplot} by number of observations with maximum conditional probability
#' greater than \code{threshold}
#' @param profiles_order If \code{NULL} or \code{FALSE}, line plots and boxplots of profiles are 
#' plotted sequentially by cluster number (K=1, K=2, ...). If \code{TRUE}, line plots and boxplots of
#' profiles are plotted in an automatically calculated order (according to the Euclidean distance
#' between cluster means) to plot clusters with similar mean profiles next to one another. 
#' Otherwise, the user may provide a vector (of length equal to the number of clusters in the 
#' given model) providing the desired order of plots. 
#' @param ...  Additional optional plotting arguments
#' 
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#' @example inst/examples/coseq-package.R
#'
#' @export
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics boxplot
#' @importFrom graphics axis
#' @importFrom graphics barplot
plot.NormMixClus <- function(x, y_profiles=NULL, K=NULL, threshold=0.8, conds=NULL,
                             average_over_conds=FALSE, 
                             graphs=c("logLike", "ICL", 
                                      "profiles", "boxplots", "probapost_boxplots",
                                      "probapost_barplots", "probapost_histogram"), 
                             order=FALSE, profiles_order=NULL, ...) {

  ## Parse ellipsis function
  arg.user <- list(...)
  if(is.null(arg.user$alpha)) arg.user$alpha<-0.3;
  
  pl_data <- data.frame(Cluster = as.numeric(substr(names(x$logLike.all), 3, 10)),
                        logLike = x$logLike.all,
                        ICL = x$ICL.all)
  ## Likelihood plot
  if("logLike" %in% graphs) {
    gg <- ggplot(pl_data, aes_string(x="Cluster", y="logLike")) + 
      geom_point() + geom_line() +
      scale_y_continuous(name = "Log-likelihood")
    print(gg)
  }
  
  ## ICL plot
  if("ICL" %in% graphs) {
    gg <- ggplot(pl_data, aes_string(x="Cluster", y="ICL")) + 
      geom_point() + geom_line() 
    print(gg)
  }
  
  ## Model-specific plots
  if("profiles" %in% graphs | "boxplots" %in% graphs | "probapost_boxplots" %in% graphs |
     "probapost_barplots" %in% graphs | "probapost_histogram" %in% graphs) {
    if("profiles" %in% graphs | "boxplots" %in% graphs) {
      if(is.null(y_profiles)) stop("y_profiles needed to plot selected graphs")
    }
    if(is.null(K) == TRUE) {
      xx <- x$ICL.results
      KK <- NULL
    }
    if(is.null(K) == FALSE) {
      if(K != "ICL") {
        if(length(which(names(x$all.results) == paste("K=",K,sep=""))) == 0)
          stop("Selected model was not estimated by coseq");
        if(length(K) > 1) 
          stop("K must be NULL, a single value, or ICL")
        xx <- x$all.results[[which(names(x$all.results) == paste("K=",K,sep=""))]]
        KK <- NULL
        if(length(xx$log.like) == 0) stop("Selected model was not estimated by coseq");
      }
      if(K == "ICL") {
        xx <- x$ICL.results;
        KK <- NULL;
      }
    }

    gr <- which(!graphs %in% c("ICL", "logLike"))
    plot(x=xx, y_profiles=y_profiles, K=KK, threshold=threshold, conds=conds,
         average_over_conds=average_over_conds, 
         graphs=graphs[gr], order = order, alpha=arg.user$alpha, 
         profiles_order=profiles_order, ...)
  }
}
