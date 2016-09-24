#' Compare corrected ICL values after data transformation
#' 
#' Compare the corrected ICL values after applying the arcsin, logit, and logMedianRef
#' transformations in a coseq analysis
#'
#' @param x A list made up of \code{coseq} objects. At the current time, this function
#' only supports the comparison of \code{coseq} objects using \code{model="Normal"} and
#' \code{transformation = c("arcsin", "logit", "logMedianRef")}
#'
#' @return A plot of corrected ICL values for the models included in \code{x} (the list
#' of \code{coseq} objects)
#' 
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#' 
#' @export
#' @example inst/examples/coseq-package.R
#'
compareICL <- function(x) {
  
  coseq_list <- x
  
  if(length(coseq_list) > 3) stop("Function only currently supported for <= 3 coseq objects")
  transf <- unlist(lapply(coseq_list, function(xx) xx$transformation))
  check <- which(transf != "arcsin" & transf != "logit" & transf != "logMedianRef")
  if(length(check) > 0) 
    stop("Function only currently supported for arcsin, logit, and logMedianRef transformations")
  ## TODO: CHECK THAT EACH TRANSFORMATION IS ONLY PRESENT ONCE
  resarcsin <- reslogit <- reslogMedianRef <- NULL
  if(length(which(transf == "arcsin")) > 0) {
    resarcsin <- coseq_list[[which(transf == "arcsin")]]
  }
  if(length(which(transf == "logit")) > 0) {
    reslogit <- coseq_list[[which(transf == "logit")]]
  }
  if(length(which(transf == "logMedianRef")) > 0) {
    reslogMedianRef <- coseq_list[[which(transf == "logMedianRef")]]
  }
  
  if(sum(resarcsin$y_profiles != reslogit$y_profiles) > 0 |
     sum(resarcsin$y_profiles != reslogMedianRef$y_profiles) > 0) 
    stop("y_profiles in coseq objects not equal -- are models estimated on same data?")
  
  ## NOTE: REPLACE 0's with smallest value > 0, 1's with largest value < 1
  
  PP <- resarcsin$y_profiles 
  PP[which(PP == 0)] <- min(PP[which(PP > 0)])
  PP[which(PP == 1)] <- max(PP[which(PP < 1)])
  
  
  n <- dim(PP)[1]
  p <- dim(PP)[2]
  qarcsin <- (n*p*log(2)) + (0.5*sum(sum(log(PP*(1-PP)))))
  qlogit <- (n*p*log(log(2))) + (sum(sum(log(PP*(1-PP)))))
  qlogmedianref <- (n*p*log(log(2))) + sum(sum(log(PP)))
  
  plotmat <- matrix(NA, nrow=0, ncol=3)
  colnames(plotmat) <- c("K", "ICL", "Transformation")
  plotmat <- rbind(plotmat, 
                   cbind(resarcsin$results$nbCluster.all, 
                         resarcsin$results$ICL.all + (2*qarcsin),
                                  rep("arcsin", length(resarcsin$results$ICL.all))),
                   cbind(reslogit$results$nbCluster.all, 
                         reslogit$results$ICL.all + (2*qlogit),
                         rep("logit", length(reslogit$results$ICL.all))),
                   cbind(reslogMedianRef$results$nbCluster.all, 
                         reslogMedianRef$results$ICL.all + (2*qlogmedianref),
                         rep("logMedianRef", length(reslogMedianRef$results$ICL.all))))
  plotdf <- data.frame(plotmat, row.names=1:nrow(plotmat), stringsAsFactors = FALSE)
  plotdf$Transformation <- factor(plotdf$Transformation)
  plotdf$ICL <-  as.numeric(plotdf$ICL)
  plotdf$K <-  as.numeric(plotdf$K)
 
  g1 <- ggplot(plotdf, aes_string(x="K", y="ICL")) + 
    geom_line(aes_string(color="Transformation")) + 
    geom_point(aes_string(color="Transformation")) +
    scale_y_continuous(name="Corrected ICL") + theme_bw()
  
  print(g1) 
}



# plot.coseq <- function(resarcsin, reslogit, reslogMedianRef, y_profiles) {
#   n=dim(y_profiles)[1]
#   p=dim(y_profiles)[2]
#   qarcsin<-(n*p*log(2)) + (0.5*sum(sum(log(y_profiles*(1-y_profiles)))))
#   qlogit<-(n*p*log(log(2))) + (sum(sum(log(y_profiles*(1-y_profiles)))))
#   qlogmedianref<-(n*p*log(log(2))) + sum(sum(log(y_profiles)))
#   
#   plot(resarcsin$nbClust,resarcsin$ICLvalue + (2*qarcsin),col="red",type="l",
#        xlab="nb cluster",ylab="ICL sur les pi")
#   points(reslogit$nbClust,reslogit$ICLvalue + (2*qlogit),col="blue",type="l")
#   points(reslogMedianRef$nbClust, reslogMedianRef$ICLvalue + (2*qlogmedianref),
#          col="magenta",type="l")
#   legend("topright",legend=c("arcsin","logit","logMedianRef"),
#          col=c("red","blue","magenta"),lty=1)
# }