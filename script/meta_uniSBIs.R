# Meta-analysis of univariate SBIs
# 2017/07/10 Crystal

setwd("~/Desktop/fishsize_varpart/")

lifehist_habitat <- read.csv("~/Desktop/fishsize_varpart/data/lifehist_habitat.csv",header = TRUE)
meanSST_F <- read.csv("~/Desktop/fishsize_varpart/output/meanSST_F.csv",header = TRUE)
cvSST_F <- read.csv("~/Desktop/fishsize_varpart/output/cvSST_F.csv",header = TRUE)

# combine the columns as variable
mean_cv <- cbind(meanSST_F[,-1],cvSST_F[,-1])
colnames(mean_cv)<-c("MeanF","MeanTemp","cvF","cvTemp")

# select the best based on total R-square
adjR2B <- read.csv("uniSBI_r2B.csv",header = TRUE)[,-1] 
adjR2B1yr <- read.csv("uniSBI_r2B_1yr.csv",header = TRUE)[,-1]
adjR2B3yr <- read.csv("uniSBI_r2B_3yr.csv",header = TRUE)[,-1]
colnames(adjR2B) <- c("mean","L95","shannon","evenness","skewness")
colnames(adjR2B1yr) <- c("mean","L95","shannon","evenness","skewness")
colnames(adjR2B3yr) <- c("mean","L95","shannon","evenness","skewness")

bestR2Bmean <- max.col(cbind(adjR2B$mean,adjR2B1yr$mean,adjR2B3yr$mean))
bestR2BL95 <- max.col(cbind(adjR2B$L95,adjR2B1yr$L95,adjR2B3yr$L95))
bestR2Bshannon <- max.col(cbind(adjR2B$shannon,adjR2B1yr$shannon,adjR2B3yr$shannon))
bestR2Bevenness <- max.col(cbind(adjR2B$evenness,adjR2B1yr$evenness,adjR2B3yr$evenness))
bestR2Bskewness <- max.col(cbind(adjR2B$skewness,adjR2B1yr$skewness,adjR2B3yr$skewness))

bestR2B <- cbind(bestR2Bmean,bestR2BL95,bestR2Bshannon,bestR2Bevenness,bestR2Bskewness)

# P value
sigF <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigF.csv",header = TRUE)[,-1] 
sigF1yr <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigF_1yr.csv",header = TRUE)[,-1] 
sigF3yr <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigF_3yr.csv",header = TRUE)[,-1] 

sigSST <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigSST.csv",header = TRUE)[,-1] 
sigSST1yr <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigSST_1yr.csv",header = TRUE)[,-1] 
sigSST3yr <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigSST_3yr.csv",header = TRUE)[,-1] 

sigSSTF <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigSSTF.csv",header = TRUE)[,-1] 
sigSSTF1yr <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigSSTF_1yr.csv",header = TRUE)[,-1] 
sigSSTF3yr <- read.csv("~/Desktop/fishsize_varpart/output/uniSBI_sigSSTF_3yr.csv",header = TRUE)[,-1] 

adjR2B_aggre <- matrix(0,nrow=27,ncol=5)

sig_Faggre <- matrix(0,nrow=27,ncol=5)
sig_SSTaggre <- matrix(0,nrow=27,ncol=5)
sig_SSTFaggre <- matrix(0,nrow=27,ncol=5)

for (i in 1:27) {
  for (j in 1:5) {
  if (bestR2B[i,j] == 1) {
  adjR2B_aggre[i,j] <- adjR2B[i,j]
  sig_Faggre[i,j] <- sigF[i,j] 
  sig_SSTaggre[i,j] <- sigSST[i,j]
  sig_SSTFaggre[i,j] <- sigSSTF[i,j]
  }
  if (bestR2B[i,j] == 2) {
    adjR2B_aggre[i,j] <- adjR2B1yr[i,j]
    sig_Faggre[i,j] <- sigF1yr[i,j] 
    sig_SSTaggre[i,j] <- sigSST1yr[i,j]
    sig_SSTFaggre[i,j] <- sigSSTF1yr[i,j]
  }
  if (bestR2B[i,j] == 3) {
    adjR2B_aggre[i,j] <- adjR2B3yr[i,j]
    sig_Faggre[i,j] <- sigF3yr[i,j] 
    sig_SSTaggre[i,j] <- sigSST3yr[i,j]
    sig_SSTFaggre[i,j] <- sigSSTF3yr[i,j]  
  }
  }
}

colnames(adjR2B_aggre) <- c("mean","L95","shannon","evenness","skewness")
colnames(sig_Faggre) <- c("mean","L95","shannon","evenness","skewness")
colnames(sig_SSTaggre) <- c("mean","L95","shannon","evenness","skewness")
colnames(sig_SSTFaggre) <- c("mean","L95","shannon","evenness","skewness")

uniSBI_sigF_aggre <- sig_Faggre
uniSBI_sigSST_aggre <- sig_SSTaggre
uniSBI_sigSSTF_aggre <- sig_SSTFaggre

write.csv(uniSBI_sigF_aggre,"~/Desktop/fishsize_varpart/output/uniSBI_sigF_aggre.csv")
write.csv(uniSBI_sigSST_aggre,"~/Desktop/fishsize_varpart/output/uniSBI_sigSST_aggre.csv")
write.csv(uniSBI_sigSSTF_aggre,"~/Desktop/fishsize_varpart/output/uniSBI_sigSSTF_aggre.csv")
