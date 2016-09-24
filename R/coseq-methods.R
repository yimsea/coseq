#' Co-expression or co-abudance analysis of high-throughput sequencing data
#' 
#' This is the primary user interface for the \code{coseq} package.
#' Generic S3 functions are implemented to perform co-expression or co-abudance analysis of 
#' high-throughput sequencing data, with or without data transformation, using mixture models. 
#' The supported classes are \code{matrix}, \code{data.frame}, \code{DESeqDataSet}, 
#' \code{DGEList}, \code{DGEExact}, \code{DGEGLM}, and \code{DGELRT}. The output of \code{coseq} 
#' is an S3 object of class \code{coseq}.
#' 
#' @inheritParams coseq_run
#' @return
#' An S3 object of class \code{coseq} containing the following:
#' \item{results }{Object of class \code{NormMixClus} or 
#' \code{PoisMixClusWrapper}, the latter being the class defined by the 
#' HTSCluster package}
#' \item{model }{Model used, either \code{Normal} or \code{Poisson}}
#' \item{transformation }{Transformation used on the data}
#' \item{tcounts }{Transformed data using to estimate model}
#' \item{y_profiles }{Normalized profiles for use in plotting}
#' 
#' @author Andrea Rau
#' @export
#' @example inst/examples/coseq-package.R
#' 
coseq <- function(y, K, subset=NULL, model="Normal", transformation="none", norm="TMM", 
                  meanFilterCutoff=NULL, modelChoice="ICL",
                  parallel=FALSE, BPPARAM=bpparam(), ...) {
  UseMethod("coseq", y)
}

## TODO check that DESeq2 and TMM norm factors match with coseq
#########################################################################################
#' @describeIn coseq Perform coseq analysis for \code{matrix} class
#' @method coseq matrix
#' @importFrom stats p.adjust
#' @export
coseq.matrix <- function(y, K, subset=NULL, model="Normal", transformation="none", 
                         norm="TMM", meanFilterCutoff=NULL, modelChoice="ICL",
                         parallel=FALSE,  BPPARAM=bpparam(), ...) {
  arg.user <- list(...)
  
  if(is.null(subset) | is.numeric(subset)) {
    run <- coseq_run(y=y, K=K, subset=subset, model=model, transformation=transformation, 
                     norm=norm, meanFilterCutoff=meanFilterCutoff, modelChoice=modelChoice,
                     parallel=parallel, 
                     BPPARAM=BPPARAM, ...)
  }
  if(class(subset) == "DESeqResults") {
    if(nrow(subset) != nrow(y)) stop("y and subset must have the same number of rows")
    res <- subset
    subset.index <- which(res$padj < res@metadata$alpha)
    cat("Co-expression analysis on DESeq2 output:\n")
    cat(paste(length(subset.index), "DE genes at p-adj <", res@metadata$alpha, "\n"))
    cat("****************************************\n")
    run <- coseq_run(y=y, K=K, subset=subset.index, model=model, transformation=transformation, 
                     norm="DESeq", meanFilterCutoff=NULL, modelChoice=modelChoice,
                     parallel=parallel, 
                     BPPARAM=BPPARAM, ...)
  }
  if(class(subset) == "DGELRT") {
    res <- subset
    if(is.null(arg.user$alpha)) {
      alpha <- 0.05
    }
    if(!is.null(arg.user$alpha)) {
      alpha <- arg.user$alpha
    }
    subset.index <- which(p.adjust(subset$table$PValue, method="BH") < alpha)
    cat("Co-expression analysis on edgeR output:\n")
    cat(paste(length(subset.index), "DE genes at p-adj <", alpha, "\n"))
    cat("****************************************\n")
    run <- coseq_run(y=y, K=K, subset=subset.index, model=model, transformation=transformation, 
                     norm="TMM", meanFilterCutoff=NULL, modelChoice=modelChoice, 
                     parallel=parallel, 
                     BPPARAM=BPPARAM, ...)
  }
  return(run)
}


#########################################################################################
#' @describeIn coseq Perform coseq analysis for \code{data.frame} class
#' @method coseq data.frame
#' @export
coseq.data.frame <- function(y, K, subset=NULL, model="Normal", transformation="arcsin", 
                             norm="TMM", meanFilterCutoff=NULL, modelChoice="ICL",
                             parallel=FALSE, BPPARAM=bpparam(), ...) {
  arg.user <- list(...)
  
  if(is.null(subset) | is.numeric(subset)) {
    run <- coseq_run(y=y, K=K, subset=subset, model=model, transformation=transformation, 
                     norm=norm, meanFilterCutoff=meanFilterCutoff, modelChoice=modelChoice,
                     parallel=parallel, 
                     BPPARAM=BPPARAM, ...)
  }
  if(class(subset) == "DESeqResults") {
    if(nrow(subset) != nrow(y)) stop("y and subset must have the same number of rows")
    res <- subset
    subset.index <- which(res$padj < res@metadata$alpha)
    cat("Co-expression analysis on DESeq2 output:\n")
    cat(paste(length(subset.index), "DE genes at p-adj <", res@metadata$alpha, "\n"))
    cat("****************************************\n")
    run <- coseq_run(y=y, K=K, subset=subset.index, model=model, transformation=transformation, 
                     norm="DESeq", meanFilterCutoff=NULL, modelChoice=modelChoice,
                     parallel=parallel, 
                     BPPARAM=BPPARAM, ...)
  }
  if(class(subset) == "DGELRT") {
    res <- subset
    if(is.null(arg.user$alpha)) {
      alpha <- 0.05
    }
    if(!is.null(arg.user$alpha)) {
      alpha <- arg.user$alpha
    }
    subset.index <- which(p.adjust(subset$table$PValue, method="BH") < alpha)
    cat("Co-expression analysis on edgeR output:\n")
    cat(paste(length(subset.index), "DE genes at p-adj <", alpha, "\n"))
    cat("****************************************\n")
    run <- coseq_run(y=y, K=K, subset=subset.index, model=model, transformation=transformation, 
                     norm="TMM", meanFilterCutoff=NULL, modelChoice=modelChoice,
                     parallel=parallel, 
                     BPPARAM=BPPARAM, ...)
  }
  return(run)
}


#########################################################################################
#' @describeIn coseq Perform coseq analysis for \code{DESeqDataSet} class from \code{DESeq2} 
#' package
#' @method coseq DESeqDataSet
#' @export
#' @importMethodsFrom DESeq2 counts sizeFactors
#' @importFrom DESeq2 results
#' 
coseq.DESeqDataSet <- function(y, K, subset=NULL, model="Normal", transformation="arcsin", 
                               norm="TMM", meanFilterCutoff=NULL, modelChoice="ICL",
                               parallel=FALSE, BPPARAM=bpparam(), ...) {
  
  ## Parse ellipsis function separately for DESeq and coseq
  dots <- dots_DESeq <- dots_coseq <- list(...)
  DESeq_argindex <- setdiff(names(dots_DESeq), c("alpha"))
  coseq_argindex <- 
    setdiff(names(dots_coseq), c("conds", "model", "transformation", "geneIDs",
                                 "parallel", "BPPARAM", "alg.type", "init.runs",
                                 "init.type", "GaussianModel", "init.iter",
                                 "cutoff", "verbose", "digits", "fixed.lambda",
                                 "equal.proportions", "prev.labels", 
                                 "prev.probaPost", "interpretation", "EM.verbose",
                                 "wrapper", "modelChoice"))
  dots_DESeq[DESeq_argindex] <- NULL
  dots_coseq[coseq_argindex] <- NULL

  count_matrix <- counts(y)
  if(is.null(subset)) {
    if(is.null(dots_DESeq$alpha)) res <- results(y);
    if(!is.null(dots_DESeq$alpha)) res <- results(y, alpha=dots_DESeq$alpha)
  }
  if(!is.null(subset)) {
    if(class(subset) != "DESeqResults") stop("subset should be the output of DESeq2 results")
    res <- subset
  }
  norm <- sizeFactors(y)
  subset.index <- which(res$padj < res@metadata$alpha)
  cat("Co-expression analysis on DESeq2 output:\n")
  cat(paste(length(subset.index), "DE genes at p-adj <", res@metadata$alpha, "\n"))
  cat("****************************************\n")
  meanFilterCutoff <- NULL
  
  run <- do.call(coseq_run, list(y=count_matrix, K=K, subset=subset.index,
                   meanFilterCutoff=meanFilterCutoff, norm=norm, dots_coseq))
  return(run)
}


# #' S4 generic method for co-expression or co-abudance analysis of high-throughput sequencing data
# #' @name coseq
# #' @rdname coseq
# #' @export
# #' @importFrom methods setGeneric setMethod
# setGeneric(
#    name = "coseq",
#    def = function(object, ...) {standardGeneric("coseq")}
# )
# 
# #' @describeIn coseq Blah blah blah
# setMethod(
#   f= "coseq",
#   signature = signature(object="matrix"),
#   definition = function(object, ...) {
#     run <- coseq_run(y=object, ...)
#     return(run)
#   }
# )