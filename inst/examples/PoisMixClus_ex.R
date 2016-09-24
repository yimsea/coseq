## Simulate toy data, n = 300 observations
set.seed(12345)
countmat <- matrix(round(runif(300*4, min=0, max=500)), nrow=300, ncol=4)
countmat <- countmat[which(rowSums(countmat) > 0),]
conds <- rep(c("A","B","C","D"), each=2)

## Run the Poisson mixture model for K = 2,3
run <- coseq(y=countmat, K=2:3, iter=5, model="Poisson")
plot(run)
summary(run)
