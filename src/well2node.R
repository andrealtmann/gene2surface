## find node that is closest to a given well_id
library(data.table)


#load lh. aibs data:
load("Data/lh_cort_fsa_newmni.RData")


rownames(sample.info.new) <- sample.info.new$well_id
use_well <- rownames(merged.dat.lh.cort)

##load node coords
target.data  <- fread("cat lh.MNI152_2mm.scannerRAS.lfname.gz | gzip -d", data.table=F)
target.coord <- target.data[,c("X","Y","Z")]

train.data  <- sample.info.new[use_well,]
train.coord <- train.data[,paste("corrected_mni",c("x","y","z"),sep="_")]

##euclidean dist
t2t <- apply(train.coord,1,function(x){
   xv <- unlist(x)
   tcdiff <- sweep(target.coord,2,xv)
   sqrt(rowSums(tcdiff^2))
})

wn.map <- t(apply(t2t,2, function(x){ idx <- which.min(x); mmm <- min(x); c(target.data[idx,"vID"],mmm) }))
