#' Permute columns of a contingency table
#' 
#' Permute the columns of a contingency table comparing two clusterings
#' to load the diagonal as much as possible.
#'
#' @param table_1 Partition from a first data clustering
#' @param table_2 Partition from a second data clustering
#'
#' @return Permuted contingency table 
#' @export
#' 
#' @examples 
#' ## Generate arbitrary labels from two separate clustering results
#' labels_1 <- sample(1:10, 1000, replace=TRUE)  ## K=10 clusters
#' labels_2 <- sample(1:8, 1000, replace=TRUE)   ## K=8 clusters
#' matchContTable(labels_1, labels_2)
#' 
#' @importFrom e1071 matchClasses
matchContTable <- function(table_1, table_2){
  tab <- table(table_1, table_2)
  ## Put larger clustering in rows if needed, nrow(tab) >= ncol(tab)
  transpose <- FALSE
  if(nrow(tab) < ncol(tab)) transpose <- TRUE;
  if(transpose==TRUE) tab <- t(tab);
  ## Order rows according to largest clusters
  tab <- tab[order(apply(tab,1,max), decreasing=TRUE),]
  ## Match best column with each row of tab
  ## Use unique indices as some columns might map to multiple rows 
  index <- matchClasses(tab, method=ifelse(nrow(tab)==ncol(tab), "exact", "rowmax"))
  tabord <- tab[,unique(index)]
  if(transpose==TRUE) tabord <- t(tabord)
  return(tabord)
}