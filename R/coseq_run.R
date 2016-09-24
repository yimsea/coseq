#' Co-expression analysis 
#'
#' Function for primary code to perform co-expression analysis, with or without data transformation,
#' using mixture models. The output of \code{coseq_run} is an S3 object of class \code{coseq}.
#' 
#' @param y (\emph{n} x \emph{q}) matrix of observed counts for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value or a vector of values)
#' @param conds Vector of length \emph{q} defining the condition (treatment
#' group) for each variable (column) in \code{y} 
#' @param norm The type of estimator to be used to normalize for differences in
#' library size: (\dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} for
#' upper quantile, \dQuote{\code{Med}} for median, \dQuote{\code{DESeq}} for
#' the normalization method in the DESeq package, and \dQuote{\code{TMM}} for
#' the TMM normalization method (Robinson and Oshlack, 2010). Can also be a
#' vector (of length \emph{q}) containing pre-estimated library size estimates
#' for each sample. 
#' @param model Type of mixture model to use (\dQuote{\code{Poisson}} or \dQuote{\code{Normal}})
#' @param transformation Transformation type to be used: \dQuote{\code{voom}}, \dQuote{\code{logRPKM}}
#' (if \code{geneLength} is provided by user), \dQuote{\code{arcsin}}, \dQuote{\code{logit}},
#' \dQuote{\code{logMedianRef}}, \dQuote{\code{profile}}, \dQuote{\code{none}}
#' @param subset Optional vector providing the indices of a subset of
#' genes that should be used for the co-expression analysis (i.e., row indices
#' of the data matrix \code{y}. For the generic function \code{coseq}, the results of a previously
#' run differential analysis may be used to select a subset of genes on which to perform the 
#' co-expression analysis. If this is desired, \code{subset.index} can also be an object of class 
#' DESeqResults (from the \code{results} function in \code{DESeq2}).
#' @param meanFilterCutoff Value used to filter low mean normalized counts if desired (by default,
#' set to a value of 50)
#' @param modelChoice Criterion used to select the best model. For Gaussian mixture models, 
#' \dQuote{\code{ICL}} (integrated completed likelihood criterion) is currently supported. For Poisson
#' mixture models, \dQuote{\code{ICL}}, \dQuote{\code{BIC}} (Bayesian information criterion), and a
#' non-asymptotic criterion calibrated via the slope heuristics  using either the \dQuote{\code{DDSE}}
#' (data-driven slope estimation) or \dQuote{\code{Djump}} (dimension jump) approaches may be used.
#' See the \code{HTSCluster} package documentation for more details about the slope heuristics approaches.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel 
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running 
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects 
#' from the current R environment before calling the function, as it is possible that R's 
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when 
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register} 
#' will be used. 
#' @param ... Additional optional parameters.
#'
#' @return
#' An S3 object of class \code{coseq} containing the following:
#' \item{results }{Object of class \code{NormMixClus} or 
#' \code{PoisMixClus}}
#' \item{model }{Model used, either \code{Normal} or \code{Poisson}}
#' \item{transformation }{Transformation used on the data}
#' \item{tcounts }{Transformed data using to estimate model}
#' \item{y_profiles }{Normalized profiles for use in plotting}
#' \item{norm }{Normalization factors used in the analysis}
#' 
#' @author Andrea Rau
#' 
#' @export
#' @importFrom HTSCluster PoisMixClus
#' @importFrom HTSCluster PoisMixClusWrapper
#' @importFrom stats na.omit
#' @importFrom capushe capushe
#'
#' @examples 
#' ## Simulate toy data, n = 300 observations
#' set.seed(12345)
#' countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
#' countmat <- countmat[which(rowSums(countmat) > 0),]
#' conds <- rep(c("A","B","C","D"), each=2)
#' 
#' ## Run the Normal mixture model for K = 2,3,4
#' ## The following are equivalent:
#' run <- coseq_run(y=countmat, K=2:4, iter=5, transformation="arcsin")
#' run <- coseq(y=countmat, K=2:4, iter=5, transformation="arcsin")
#'
#' 
coseq_run <- function(y, K, conds=NULL, norm="TMM", model="Normal", transformation="arcsin", 
                       subset=NULL, meanFilterCutoff=50, modelChoice="ICL",
                       parallel=FALSE, BPPARAM=bpparam(), ...) {
  
  subset.index <- subset
  
  ## Parse ellipsis function
  arg.user <- list(...)

  ## Optional parameters for PoisMixClus
  if(model == "Poisson") {
    if(is.null(arg.user$Kmin.init)) arg.user$Kmin.init<-"small-em";
    if(is.null(arg.user$split.init)) arg.user$split.init<-FALSE;
    if(is.null(arg.user$init.runs)) arg.user$init.runs<-1;
    if(is.null(arg.user$init.iter)) arg.user$init.iter<-10;
    if(is.null(arg.user$alg.type)) arg.user$alg.type<-"EM";
    if(is.null(arg.user$cutoff)) arg.user$cutoff<-1e-05;
    if(is.null(arg.user$iter)) arg.user$iter<-1000;
    if(is.null(arg.user$fixed.lambda)) arg.user$fixed.lambda<-NA;
    if(is.null(arg.user$equal.proportions)) arg.user$equal.proportions<-FALSE;
    if(is.null(arg.user$verbose)) arg.user$verbose<-FALSE;
    if(is.null(arg.user$EM.verbose)) arg.user$EM.verbose<-FALSE;
    if(is.null(arg.user$interpretation)) arg.user$interpretation<-"sum";
    if(is.null(arg.user$digits)) arg.user$digits<-3; 
    if(is.null(subset)) subset.index <- NA
  }

  ## Optional parameters for transform_RNAseq
  if(is.null(arg.user$geneLength)) arg.user$geneLength<-NA;
  
  ## Optional parameters for NormMixClus
  if(model == "Normal") {
    if(is.null(arg.user$alg.type)) arg.user$alg.type<-"EM";
    if(is.null(arg.user$init.runs)) arg.user$init.runs<-50;
    if(is.null(arg.user$init.type)) arg.user$init.type<-"small-em";
    if(is.null(arg.user$init.iter)) arg.user$init.iter<-20;
    if(is.null(arg.user$iter)) arg.user$iter<-1000;
    if(is.null(arg.user$cutoff)) arg.user$cutoff<-0.001;
    if(is.null(arg.user$digits)) arg.user$digits<-3; 
  }
  
  y_profiles <- round(transform_RNAseq(y=y, norm=norm, transformation="profile",
                                       meanFilterCutoff=meanFilterCutoff, verbose=TRUE)$tcounts,
                      digits=arg.user$digits)
  
  ########################
  ## POISSON MIXTURE MODEL
  ########################
  if(length(model) == 1 & model == "Poisson") {
    if(transformation != "none") stop("Poisson mixture model may only be applied on raw counts.")
    if(is.null(conds)) {
      message("Poisson mixture model fit assuming each sample is an independent condition.")
      conds <- 1:ncol(y)
    }
    
    ## Grouping columns of y in order of condition (all replicates put together)
    o.ycols <- order(conds)
    y <- y[,o.ycols]
    conds <- conds[o.ycols]
    conds.names <- unique(conds)
    d <- length(unique(conds))
    r <- as.vector(table(conds))
    if(length(rownames(y)) == 0) rn <- 1:nrow(y);
    if(length(rownames(y)) > 0) rn <- rownames(y);
    y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
    rownames(y) <- rn;
    
    ## In case only a subset of data are to be used for analysis
    if(!is.null(subset)) {
      y <- y[subset.index,]
      n <- dim(y)[1];cols <- dim(y)[2]
      w <- rowSums(y)
    }
    if(is.null(subset)) {
      n <- dim(y)[1];cols <- dim(y)[2]
      w <- rowSums(y)
    }
    
    tcounts <- transform_RNAseq(y=y, norm=norm, transformation="none", 
                                geneLength=arg.user$geneLength, 
                                meanFilterCutoff=meanFilterCutoff, verbose=FALSE)
    
    if(parallel == FALSE) {
      run <- suppressWarnings(PoisMixClusWrapper(y=tcounts$tcounts, gmin=min(K), 
                                                 gmax=max(K), conds=conds,
                                norm=tcounts$ellnorm / sum(tcounts$ellnorm), 
                                gmin.init.type=arg.user$Kmin.init, 
                                split.init=arg.user$split.init, subset.index=subset.index, 
                                init.runs=arg.user$init.runs, init.iter=arg.user$init.iter, 
                                alg.type=arg.user$alg.type, cutoff=arg.user$cutoff, 
                                iter=arg.user$iter,
                                fixed.lambda=arg.user$fixed.lambda, 
                                equal.proportions=arg.user$equal.proportions,
                                verbose=arg.user$verbose, EM.verbose=arg.user$EM.verbose, 
                                interpretation=arg.user$interpretation))
      names(run$all.results) <- paste("K=", K, sep="")
      run$nbCluster.all <- K
    }
    
    
    if(parallel == TRUE) {
      all.results <- vector("list", length = length(K))
      names(all.results) <- paste("K=", K, sep = "")
      cat("Running K =", min(K), "...\n")
      if(arg.user$split.init == TRUE) {
        warning("Splitting initialization is not compatible with parallelization.")
      }
      run <- PoisMixClus(y=tcounts$tcounts, g=min(K), conds=conds, 
                         norm=tcounts$ellnorm / sum(tcounts$ellnorm), 
                         init.type=arg.user$Kmin.init,  
                         subset.index=subset.index, 
                         wrapper=TRUE, init.runs=arg.user$init.runs,
                         init.iter=arg.user$init.iter, alg.type=arg.user$alg.type,
                         cutoff=arg.user$cutoff, iter=arg.user$iter, 
                         fixed.lambda=arg.user$fixed.lambda,
                         equal.proportions=arg.user$equal.proportions,
                         prev.labels=NA, prev.probaPost = NA,
                         verbose=arg.user$verbose,
                         EM.verbose=arg.user$EM.verbose, 
                         interpretation=arg.user$interpretation)
      all.results[[1]] <- run
      
      index <- 2
      remainingK <- K[-which(K == min(K))]
      if(length(remainingK) > 0) {
        tmp <- bplapply(remainingK, function(ii) {
          cat("Running K =", ii, "...\n")
          res <- PoisMixClus(g=as.numeric(ii), y=tcounts$tcounts, 
                             conds=conds, norm=tcounts$ellnorm / sum(tcounts$ellnorm), 
                             init.type=arg.user$Kmin.init, 
                             subset.index=subset.index, wrapper=TRUE, 
                             prev.probaPost=NA, prev.labels=NA, 
                             init.runs = arg.user$init.runs, init.iter = arg.user$init.iter, 
                             alg.type = arg.user$alg.type, cutoff = arg.user$cutoff,
                             iter = arg.user$iter, fixed.lambda = arg.user$fixed.lambda, 
                             equal.proportions = arg.user$equal.proportions,
                             verbose = arg.user$verbose,
                             interpretation = arg.user$interpretation, 
                             EM.verbose = arg.user$EM.verbose)
          return(res)}, BPPARAM=BPPARAM)
        Kmods <- paste("K=", unlist(lapply(tmp, function(x) ncol(x$lambda))), sep="")
        all.results[-1] <- tmp[na.omit(match(names(all.results), Kmods))]
      }
      
      logLike.all <- unlist(lapply(all.results, function(x) x$log.like))
      ICL.all <- unlist(lapply(all.results, function(x) x$ICL))
      ICL.choose <- which(ICL.all == max(ICL.all, na.rm = TRUE))
      select.results <- all.results[[ICL.choose]]
      select.results$model.selection <- "ICL"
      
      BIC.all <- unlist(lapply(all.results, function(x) x$BIC))
      BIC.choose <- which(BIC.all == max(BIC.all, na.rm = TRUE))
      select.results2 <- all.results[[BIC.choose]]
      select.results2$model.selection <- "BIC"
      
      # Apply capushe: only if at least 10 models are considered
      if(c(max(K) - min(K) + 1) <= 10) {
        message("Note: slope heuristics for model selection only applied if > 10 models are fit.")
        DDSE.results <- NA
        Djump.results <- NA
        capushe <- NA
        ResCapushe <- NA
      }
      if(c(max(K) - min(K) + 1) > 10) {
        message("Note: diagnostic plots for slope heuristics (Djump and DDSE) should be examined to ensure that sufficiently complex models have been considered.")
        Kchoice <- K
        np <- (Kchoice-1) + (length(unique(conds))-1)*(Kchoice)
        mat <- cbind(Kchoice, np/n, np/n, -logLike.all/n)
        ResCapushe <- suppressWarnings(capushe(mat, n))
        DDSE <- ResCapushe@DDSE@model
        Djump <- ResCapushe@Djump@model
        DDSE.results <- all.results[[paste("K=", DDSE, sep="")]]
        Djump.results <- all.results[[paste("K=", Djump, sep="")]]
        DDSE.results$model.selection <- "DDSE"
        Djump.results$model.selection <- "Djump"
      }
      
      run <- list(nbCluster.all=K, logLike.all = logLike.all, 
                  ICL.all = ICL.all,
                  capushe = ResCapushe, 
                  all.results = all.results,
                  DDSE.results = DDSE.results,
                  Djump.results = Djump.results,
                  BIC.results = select.results2,
                  ICL.results = select.results)
      class(run) <- "HTSClusterWrapper"
    }
    
    if(modelChoice == "BIC") final.results <- run$BIC.results
    if(modelChoice == "ICL") final.results <- run$ICL.results
    if(modelChoice == "DDSE") final.results <- run$DDSE.results
    if(modelChoice == "Djump") final.results <- run$Djump.results
    
    final.results$probaPost <- round(final.results$probaPost, arg.user$digits)
    for(jj in 1:length(run$all.results)) {
      run$all.results[[jj]]$probaPost <- round(run$all.results[[jj]]$probaPost,
                                               arg.user$digits)
    }
    
    final.run <- list(nbCluster.all=run$nbCluster.all, logLike.all=run$logLike.all,
                      ICL.all=run$ICL.all, selected.results=final.results, 
                      all.results=run$all.results,
                      capushe=run$capushe)
    run <- final.run
    class(run) <- "PoisMixClus"
  }
  
  
  
  ########################
  ## NORMAL MIXTURE MODEL
  ########################
  if(length(model) == 1 & model == "Normal") {
    
    tcounts <- transform_RNAseq(y=y, norm=norm, transformation=transformation, 
                                geneLength=arg.user$geneLength,
                                meanFilterCutoff=meanFilterCutoff, verbose=FALSE)
#    if(parallel == TRUE) arg.user$verbose <- FALSE;
    run <- NormMixClus(y_profiles=tcounts$tcounts, K=K, subset=subset.index, 
                       parallel=parallel,
                       BPPARAM=BPPARAM, alg.type=arg.user$alg.type, 
                       init.runs=arg.user$init.runs,
                       init.type=arg.user$init.type, init.iter=arg.user$init.iter, 
                       iter=arg.user$iter, cutoff=arg.user$cutoff,
                       verbose=arg.user$verbose, digits=arg.user$digits)
  }

  ####################################
  ## POISSON AND NORMAL MIXTURE MODELS? To be added later...
  ####################################
  # if("Normal" %in% model & "Poisson" %in% model) {
  #    
  # }
  
  
  ####################################
  ## RETURN RESULTS
  ####################################
  if(!is.null(subset)) {
    tcounts$tcounts <- tcounts$tcounts[subset.index,]
    y_profiles <- y_profiles[subset.index,]
  }
  
  RESULTS <- list(results=run, model=model, transformation=transformation, 
                  tcounts=tcounts$tcounts, y_profiles=y_profiles, norm=tcounts$snorm)
  class(RESULTS) <- "coseq"
  return(RESULTS)
}