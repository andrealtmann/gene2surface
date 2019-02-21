## creates predictions for every node
## for a specified probe_id
library(data.table)

source("src/helper_functions.R")

##settings

probe_id <- "A_23_P48325" #KCTD4 

#uses corrected MNI coorinates from: https://github.com/chrisfilo/alleninf
dist_pref <- "corrected_mni"

#distance type
#dtype <- "eucl"
dtype <- "geod"

ofname <- paste("results/mean_", dtype, "_", probe_id, ".txt", sep="")

sigma.range <- seq(1,10,0.25)
nboot <- 100

##
exp.fname <- "Data/lh_cort_fsa_newmni.RData"
geo.fname <- "Data/mni_2mm_geo_newmni.RData"
tar.fname <- "Data/lh.MNI152_2mm.scannerRAS.lfname.gz"

if (!exists("geo.loaded"))
  geo.loaded <- F

message("load expression data")
load(exp.fname)
rownames(sample.info.new) <- sample.info.new$well_id
sample.info.new <- sample.info.new[rownames(merged.dat.lh.cort),]
train.data <- data.frame(sample.info.new, gene=merged.dat.lh.cort[,probe_id])

message("load target coordinates")
target.data  <- fread(paste("cat ",tar.fname," | gzip -d"), data.table=F)
target.coord <- target.data[,c("X","Y","Z")]
idx.node <- target.data[,"vID"] + 1

if (dtype == "eucl"){
  message("computing euclidean distances")
  train.coord <- train.data[,paste(dist_pref,c("x","y","z"),sep="_")]

  train2train <- as.matrix(dist(train.coord))
  train2test  <- apply(train.coord, 1, function(x){
    xv <- unlist(x)
    tcdiff <- sweep(target.coord, 2, xv)
    sqrt(rowSums(tcdiff^2))
  })
}

if (dtype == "geod" & !geo.loaded){
  message("load geo distances")
  load(geo.fname)
  geo.loaded <- T
  worder <- rownames(train.data)

  rownames(w2n) <- w2n$well_id
  w2n  <- w2n[worder,]
  geod <- geod[worder,]

  ##add eucledian dist from sample to closest node into the mix
  train2test <- t(sweep(geod,1, w2n[,"edist"], FUN="+"))
  #train2test <- t(sweep(geod,1, w2n[,"edist"], FUN="+"))

  nix <- w2n$node_id + 1
  train2train <- as.matrix(train2test[nix,])
  diag(train2train) <- 0

  train2test <- train2test[idx.node,]

}


message("computing optimal sigma for each donor")
donors <- paste(unique(train.data$sampleID))

donor.sig <- sapply(donors, function(x){

  idx <- train.data$sampleID==x
  tmp.t2t <- train2train[idx,idx]
  tmp.exp <- train.data[idx,"gene"]

  message("running bootstrap estimate for donor: ", x)
  boot <- bootstrap.pred(tmp.t2t, tmp.exp, sigma.range, nboot)
  tmp <- apply(boot,2,mean)
  sig.id <- which.min(tmp)
  return(sigma.range[sig.id])

})

message("making predictions")
#now make predictions for each vertex!
donor.pred <- sapply(donors, function(x){
  idx <- train.data$sampleID==x
  tmp.tr2te <- dnorm(train2test[,idx],sd=donor.sig[x])
  tmp.exp <- train.data[idx,"gene"]

  pv <- apply(tmp.tr2te,1,function(w){
    weighted.mean(tmp.exp,w)
  })
})

##prepare output
donor.mean <- apply(donor.pred, 1, mean)

max.node <- max(idx.node)

to.write <- rep(0, max.node)
to.write[idx.node] <- donor.mean
write.table(to.write, ofname, row.names=F, col.names=F, quote=F)
