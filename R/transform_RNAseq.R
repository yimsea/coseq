#' Transform RNA-seq data using common transformations
#' 
#' Application of common transformations for RNA-seq data prior to fitting a normal mixture model
#'
#' @param y (\emph{n} x \emph{q}) matrix of observed counts for \emph{n}
#' observations and \emph{q} variables
#' @param norm The type of estimator to be used to normalize for differences in
#' library size: \dQuote{\code{TC}} for total count, \dQuote{\code{DESeq}} for
#' the normalization method in the DESeq package, and \dQuote{\code{TMM}} for
#' the TMM normalization method (Robinson and Oshlack, 2010). Can also be a
#' vector (of length \emph{q}) containing pre-estimated library size estimates
#' for each sample. 
#' @param transformation Transformation type to be used: \dQuote{\code{arcsin}}, 
#' \dQuote{\code{logit}}, \dQuote{\code{logMedianRef}}, \dQuote{\code{profile}}, 
#' \dQuote{\code{voom}}, \dQuote{\code{logRPKM}} (if \code{geneLength} is provided by user),
#' \dQuote{\code{none}}
#' @param geneLength Vector of length equal to the number of rows in \dQuote{\code{y}} providing
#' the gene length (bp) for RPKM calculation
#' @param meanFilterCutoff Value used to filter low mean normalized counts 
#' @param verbose If \code{TRUE}, include verbose output
#'
#' @return
#' \item{tcounts }{Transformed counts}
#' \item{normCounts }{Normalized counts}
#' \item{snorm }{Per-sample normalization factors divided by mean normalization factor}
#' \item{ellnorm }{Per-sample normalization factors}
#' 
#' @export
#' 
#' @examples 
#' set.seed(12345)
#' countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
#' countmat <- countmat[which(rowSums(countmat) > 0),]
#' conds <- rep(c("A","B","C","D"), each=2)
#' 
#' ## Arcsin transformation, TMM normalization
#' arcsin <- transform_RNAseq(countmat, norm="TMM", transformation="arcsin")$tcounts
#' ## Logit transformation, TMM normalization
#' logit <- transform_RNAseq(countmat, norm="TMM", transformation="logit")$tcounts
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR cpm
#' @importFrom HTSFilter HTSBasicFilter
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom stats median


transform_RNAseq <- function(y, norm="TMM", transformation="arcsin", 
                             geneLength=NA, meanFilterCutoff=NULL, verbose=TRUE) {

  ##################################
  ## Calculate normalization factors
  ##################################
  ## Only calculate s values if they are not provided
  if(length(norm) != 1 & length(norm) != ncol(y)) stop(paste(sQuote("norm"), "must be one of
           the following: none, TC, DESeq, TMM, or a vector oflength equal to the number of columns 
                                                                   in", sQuote("y")))
  ## If estimated from data, all genes should be used
  if(length(norm) == 1) {
    if(norm == "none") {
      libsize <- rep(1, ncol(y))
    }
    if(norm == "TMM") {
      f <- calcNormFactors(as.matrix(y), method="TMM")
      libsize <- colSums(y) * f
    }
    if(norm == "TC") {
      libsize <- colSums(y)
    }
    if(norm == "DESeq") {
      loggeomeans <- rowMeans(log(y))
      f <- apply(y, 2, function(x) 
        exp(median((log(x) - loggeomeans)[is.finite(loggeomeans)])))
      libsize <- colSums(y) * f
    } 
  }
  if(length(norm) > 1) {
    libsize <- norm
  }

  snorm <- libsize / (mean(libsize)) 
  ellnorm <- libsize
  
  ##################################
  ## Filter data based on mean count if desired
  ##################################
  if(is.null(meanFilterCutoff) == FALSE) {

      filter <- HTSBasicFilter(y, method="mean", cutoff.type="value", 
                               cutoff=meanFilterCutoff, normalization=norm)
      if(verbose == TRUE) {
        cat("Filter applied: remove observations with normalized mean < ", meanFilterCutoff, "\n")
        cat("               ", nrow(filter$filteredData), "observations retained for analysis\n")
      }
      y <- filter$filteredData
      geneLength <- unlist(geneLength)[which(filter$on == 1)]
  }
  
  ##################################
  ## Transform data and calculate the first derivative of transformations
  ##################################
  
  normCounts <- t(t(y)/snorm +1)
  if(transformation == "none") {
    tcounts <- y
  }
  if(transformation == "profile") {
    tcounts <- normCounts / rowSums(normCounts)
  }
  if(transformation == "clrProfile") {
    profiles <- normCounts / rowSums(normCounts)
    gm <- apply(profiles, 1, function(x) exp(mean(log(x))))
    tcounts <- log(profiles / gm)
  }
  if(transformation == "voom") {
    tcounts <- log2(t((t(y + 0.5)/(ellnorm + 1)))*1e+06)
  }
  if(transformation == "logRPKM") {
    if(is.na(geneLength)[1]==TRUE) stop("RPKM transformation requires input for gene lengths")
    if(length(unlist(geneLength))!=nrow(y)) stop("Gene length vector of different dimension than counts")
    Lcounts <- y/geneLength
    tcounts <- log2(t(t(Lcounts)/ellnorm) * 10^9 + 0.5)
    elltmp <- matrix(rep(ellnorm, each = nrow(y)), nrow=nrow(y))
    lentmp <- matrix(rep(geneLength, times=ncol(y)), nrow=nrow(y))
  }
  if(transformation == "arcsin") {
    props <- normCounts / rowSums(normCounts)
    tcounts <- asin(sqrt(props))
    w <- matrix(rep(rowSums(normCounts), times=ncol(y)), nrow=nrow(y))
  }
  if(transformation == "logit") {
    props <- normCounts / rowSums(normCounts)
    tcounts <- log2(props / (1-props))
    w <- matrix(rep(rowSums(normCounts), times=ncol(y)), nrow=nrow(y))
    stemp <- matrix(rep(snorm, each = nrow(y)), nrow=nrow(y))
  }
  if(transformation == "logMedianRef") {
    m <- apply((t(t(y)/snorm))+1, 1, median)
    tcounts <- log2((t(t(y)/snorm) + 1)/(m+1))
    stemp <- matrix(rep(snorm, each = nrow(y)), nrow=nrow(y))
  }
  ##################################
  ## Old transformations (kept for reference)
  ##################################
  if(transformation == "log") { 
    tcounts <- log2(y + 1) 
  }
  if(transformation == "normlog") {
    tcounts <- log2(t(t(y)/snorm) + 1)  
  }
  if(transformation == "logMeanRef") {
    m <- apply((t(t(y)/snorm))+1, 1, mean)
    tcounts <- log2((t(t(y)/snorm) + 1)/(m+1))
  }
  if(transformation == "logGMeanRef") {
    m <- apply((t(t(y)/snorm))+1, 1, function(x) (prod(x+1))^(1/length(x)))
    tcounts <- log2((t(t(y)/snorm) + 1)/(m))
  }
  if(transformation == "vst") {
    tcounts <- varianceStabilizingTransformation(as.matrix(y), blind=TRUE, 
                                                 fitType="parametric")
  }
  if(transformation == "moderatedCPM") {
    tcounts <- cpm(as.matrix(y), normalized.lib.sizes=TRUE, log=TRUE,
                   prior.count=0.25)  
  }
  return(list(tcounts=tcounts, normCounts=normCounts, snorm=snorm, ellnorm=ellnorm))
}

