#' Visualize results from coseq clustering
#' 
#' Plot a coseq object.
#'
#' @param x An object of class \code{"coseq"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables to be used for graphing results (optional for
#' \code{logLike}, \code{ICL}, \code{probapost_boxplots}, and \code{probapost_barplots},
#' and by default takes value \code{x$tcounts} if \code{NULL})
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
#' conditional probabilities over all clusters), \code{"lambda_barplots"} (barplots of
#' estimated proportions of counts per condition in each cluster for the Poisson mixture
#' model)
#' @param order If \code{TRUE}, order clusters in \code{probapost_boxplot} by median and
#' \code{probapost_barplot} by number of observations with maximum conditional probability
#' greater than \code{threshold}
#' @param profiles_order If \code{NULL} or \code{FALSE}, line plots and boxplots of profiles are 
#' plotted sequentially by cluster number (K=1, K=2, ...). If \code{TRUE}, line plots and boxplots of
#' profiles are plotted in an automatically calculated order (according to the Euclidean distance
#' between cluster means) to plot clusters with similar mean profiles next to one another. 
#' Otherwise, the user may provide a vector (of length equal to the number of clusters in the 
#' given model) providing the desired order of plots.
#' @param n_row Number of rows for plotting layout of line plots and boxplots of profiles.
#' Note that if \code{n_row} x \code{n_col} is less than the total number of clusters plotted,
#' plots will be divided over multiple pages.
#' @param n_col Number of columns for plotting layout of line plots and boxplots of profiles.
#' Note that if \code{n_row} x \code{n_col} is less than the total number of clusters plotted,
#' plots will be divided over multiple pages.
#' @param ...  Additional optional plotting arguments
#' 
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#' @example inst/examples/coseq-package.R
#'
#' @export
#' 
plot.coseq <- function(x, y_profiles=NULL, K=NULL, threshold=0.8, conds=NULL,
                             average_over_conds=FALSE, 
                             graphs=c("logLike", "ICL", 
                                      "profiles", "boxplots", "probapost_boxplots",
                                      "probapost_barplots", "probapost_histogram",
                                      "lambda_barplots"), 
                             order=FALSE, profiles_order=NULL, n_row=NULL, n_col=NULL, ...) {
  
  if(is.null(y_profiles) == TRUE) y_profiles <- x$y_profiles

  plot(x$results, y_profiles=y_profiles, K=K, threshold=threshold, conds=conds,
       average_over_conds=average_over_conds, graphs=graphs, 
       order=order, profiles_order=profiles_order, n_row=n_row, n_col=n_col, ...)
  
}
