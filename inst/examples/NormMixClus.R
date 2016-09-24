## Simulate toy data, n = 300 observations
set.seed(12345)
countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
countmat <- countmat[which(rowSums(countmat) > 0),]
profiles <- transform_RNAseq(countmat, norm="none", 
                             transformation="arcsin")$tcounts

conds <- rep(c("A","B","C","D"), each=2)

## Run the Normal mixture model for K = 2,3
run <- NormMixClus(y=profiles, K=2:3, iter=5)

## Run the Normal mixture model for K=2
run2 <- NormMixClus_K(y=profiles, K=2, iter=5)

## Re-estimate mixture parameters for the model with K=2 clusters
param <- NormMixParam(run2, y_profiles=profiles)

## Summary of results
summary(run, y_profiles=profiles)
