#some functions

weighted.mean <- function(vals, w){
  res <- sum(vals * w)/sum(w)
  return(res)
}

##bootstrap predictions
bootstrap.pred <- function(mdist, mval, smooths, nboot=100, frac=0.8){

  res <- c()
  N <- length(mval)
  Ntrain <- floor(N*frac)
  for(b in 1:nboot){
    #sample
    idx.train <- sample(1:N, Ntrain)
    idx.test  <- setdiff(1:N, idx.train)


    tmp <- sapply(smooths, function(si){
      tmpdist <- dnorm(mdist, sd=si)
      pv <- sapply(idx.test, function(i){
         w  <- tmpdist[i,idx.train]
         return(weighted.mean(mval[idx.train],w))
      })
      #compute error
      rmse <- sqrt(mean((pv - mval[idx.test])^2))
      return(rmse)
    })
    res <- rbind(res, tmp)
  }
  return(res)

}
