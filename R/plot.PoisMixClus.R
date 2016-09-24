#' Visualize results from clustering using a Poisson mixture model
#' 
#' Plot a PoisMixClus object.
#'
#' @param x An object of class \code{"PoisMixClus"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables to be used for graphing results (optional for
#' \code{logLike}, \code{ICL}, \code{probapost_boxplots}, and \code{probapost_barplots})
#' @param K If desired, the specific model to use for plotting. If \code{NULL},
#' the model chosen using the criterion defined by the user will be plotted
#' @param threshold Threshold used for maximum conditional probability; only observations
#' with maximum conditional probability greater than this threshold are visualized 
#' @param conds Condition labels
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
#' @param ...  Additional optional plotting arguments
#' 
#' @author Andrea Rau
#'
#' @export
#' @example inst/examples/PoisMixClus_ex.R
#' @importFrom graphics plot points boxplot axis barplot legend par mtext text layout abline
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats as.dendrogram
#' 
plot.PoisMixClus <- function(x, y_profiles=NULL, K=NULL, threshold=0.8, conds=NULL,
                             average_over_conds=FALSE, 
                             graphs=c("logLike", "ICL", 
                                      "profiles", "boxplots", "probapost_boxplots",
                                      "probapost_barplots", "probapost_histogram",
                                      "lambda_barplots"), 
                             order=FALSE, profiles_order=NULL, ...) {

  ## Parse ellipsis function
  arg.user <- list(...)
  if(is.null(arg.user$alpha)) arg.user$alpha<-0.3;
  
  pl_data <- data.frame(Cluster = as.numeric(substr(names(x$logLike.all), 3, 10)),
                        logLike = x$logLike.all,
                        ICL = -x$ICL.all)
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
     "probapost_barplots" %in% graphs | "probapost_histogram" %in% graphs | 
     "lambda_barplots" %in% graphs) {
    if("profiles" %in% graphs | "boxplots" %in% graphs) {
      if(is.null(y_profiles)) stop("y_profiles needed to plot selected graphs")
    }
    if(is.null(K) == TRUE) {
      xx <- x$selected.results
      KK <- NULL
    }
    if(is.null(K) == FALSE) {
      if(length(which(names(x$all.results) == paste("K=",K,sep=""))) == 0)
        stop("Selected model was not estimated by coseq");
      if(length(K) > 1) 
        stop("K must be NULL or a single value")
      xx <- x$all.results[[which(names(x$all.results) == paste("K=",K,sep=""))]]
      KK <- NULL
      if(length(xx$log.like) == 0) stop("Selected model was not estimated by coseq");
    }
    
    ## Change class to use NormMixClus_K plotting functions
    if("profiles" %in% graphs | "boxplots" %in% graphs | "probapost_boxplots" %in% graphs |
      "probapost_barplots" %in% graphs | "probapost_histogram" %in% graphs) {
      xx_nmc <- xx
      class(xx_nmc) <- "NormMixClus_K"
      gr <- which(!graphs %in% c("ICL", "logLike", "lambda_barplots"))
      plot(x=xx_nmc, y_profiles=y_profiles, K=KK, threshold=threshold, conds=conds,
           average_over_conds=average_over_conds, 
           graphs=graphs[gr], order = order, alpha=arg.user$alpha, 
           profiles_order=profiles_order, ...) 
    }
    
    ## Add lambda barplots if desired
    if("lambda_barplots" %in% graphs) {
 
     lambda <- xx$lambda
      pi <- xx$pi
      s <- xx$norm
      conds <- xx$conds
      ss <- rep(NA, length(unique(conds)))
      for (c in 1:length(unique(conds))) {
        ss[c] <- sum(s[which(conds == unique(conds)[c])])
      }
      lams <- lambda*ss

      ## Hclust on distance between s*lambda to order bars
      lamshc <- hclust(dist(t(lams)))
      lamsden <- as.dendrogram(lamshc)
      lmat <- cbind(rbind(c(NA, 3), 2:1, c(NA, 4)), c(0,5,0))
      lwid <- c(0.1,3.75,1)
      lhei <- c(1+0.2, 3.5, 0.5)
      lmat[is.na(lmat)] <- 0
      
      cols <- brewer.pal(max(3, length(unique(conds))),"Greens")
      if(length(unique(conds)) == 2) cols <- cols[c(1,3)]
      names(cols) <- unique(conds)
      
      layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
      ## Barplot
      par(mar=rep(0,4))
      b <- barplot(lams[,lamshc$order], beside = FALSE, col = cols[1:length(unique(conds))],
                   yaxt="n", xaxt="n", space=0)
      if (length(unique(conds)) == 2) abline(h = 0.5, lty = 2, col = "grey");
      ## Y axis label
      par(mar=c(0,0,0,0))
      plot(0,1, xaxt="n", col="white", yaxt="n", ylim=c(0,1), xlab="", frame.plot="F")
      axis(side=2, at=seq(-0.03,1.03,length=6), labels=seq(0,1,by=0.2))
      mtext(side = 2, outer = FALSE, expression(paste(lambda[jk], s[j.])), line = 2)
      ## Dendrogram
      par(mar=c(0,1,0,1))
      plot(lamsden , axes = FALSE, xaxs = "i", leaflab = "none")
      ## Cluster labels
      par(mar=c(1,1,0,1))
      plot(1:length(lamshc$order), rep(.8, length(lamshc$order)), col="white", ylim=c(0,1), 
           bty="n",
           xaxt="n", yaxt="n")
      text(1:length(lamshc$order), rep(.8, length(lamshc$order)), label = lamshc$order, 
           cex=1.25)
      points(1:length(lamshc$order), rep(.3, length(lamshc$order)), cex=10*pi, col="black", 
             pch = 16)
      mtext(side=1, "Cluster number and relative size", line = 0)
      ## Sample labels
      par(mar=c(0,1,0,0))
      plot(0, xaxt="n", yaxt="n", bty="n", col="white")
      legend(rownames(lambda), x="center", bty="n", fill=cols[1:length(unique(conds))], cex=1, 
             x.intersp=0.25)
    }
  }
  # if("capushe_diagnostics" %in% graphs) {
  #   if(!is.null(x$capushe)) {
  #     plot(x$capushe, ask=FALSE)
  #   }
  # }
}
