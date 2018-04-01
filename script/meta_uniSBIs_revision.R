# Meta-analysis of univariate SBIs
# 2017/07/10 Crystal
# 2018/02/28 Revise with new survey data analysis result
lifehist_habitat <- read.csv("~/Desktop/fishsize_varpart/data/lifehist_habitat.csv",header = TRUE)
meanSST_F <- read.csv("output/meanSST_F.csv",header = TRUE)
cvSST_F <- read.csv("output/cvSST_F.csv",header = TRUE)

# combine the columns as variable
mean_cv <- cbind(meanSST_F[,-1],cvSST_F[,-1])
colnames(mean_cv)<-c("MeanF","MeanTemp","cvF","cvTemp")

# select the best based on total R-square
adjR2B <- read.csv("output/uniSBI_r2B.csv",header = TRUE)[,-1] 
adjR2B[adjR2B<0] <-0
adjR2B1yr <- read.csv("output/uniSBI_r2B1yr.csv",header = TRUE)[,-1]
adjR2B1yr[adjR2B1yr<0] <-0
adjR2B3yr <- read.csv("output/uniSBI_r2B3yr.csv",header = TRUE)[,-1]
adjR2B3yr[is.na(adjR2B3yr)] <- 0 #replace NA with 0
adjR2B3yr[adjR2B3yr<0] <-0

colnames(adjR2B) <- c("mean","L95","shannon","evenness","skewness")
colnames(adjR2B1yr) <- c("mean","L95","shannon","evenness","skewness")
colnames(adjR2B3yr) <- c("mean","L95","shannon","evenness","skewness")

bestR2Bmean <- max.col(cbind(adjR2B$mean,adjR2B1yr$mean,adjR2B3yr$mean), ties.method = "first")
bestR2BL95 <- max.col(cbind(adjR2B$L95,adjR2B1yr$L95,adjR2B3yr$L95), ties.method = "first")
bestR2Bshannon <- max.col(cbind(adjR2B$shannon,adjR2B1yr$shannon,adjR2B3yr$shannon), ties.method = "first")
bestR2Bevenness <- max.col(cbind(adjR2B$evenness,adjR2B1yr$evenness,adjR2B3yr$evenness), ties.method = "first")
bestR2Bskewness <- max.col(cbind(adjR2B$skewness,adjR2B1yr$skewness,adjR2B3yr$skewness), ties.method = "first")

bestR2B <- cbind(bestR2Bmean,bestR2BL95,bestR2Bshannon,bestR2Bevenness,bestR2Bskewness)

# P value
sigF <- read.csv("output/uniSBI_sigF.csv",header = TRUE)[,-1] 
sigF1yr <- read.csv("output/uniSBI_sigF1yr.csv",header = TRUE)[,-1] 
sigF3yr <- read.csv("output/uniSBI_sigF3yr.csv",header = TRUE)[,-1] 

sigSST <- read.csv("output/uniSBI_sigSST.csv",header = TRUE)[,-1] 
sigSST1yr <- read.csv("output/uniSBI_sigSST1yr.csv",header = TRUE)[,-1] 
sigSST3yr <- read.csv("output/uniSBI_sigSST3yr.csv",header = TRUE)[,-1] 

sigSSTF <- read.csv("output/uniSBI_sigSSTF.csv",header = TRUE)[,-1] 
sigSSTF1yr <- read.csv("output/uniSBI_sigSSTF1yr.csv",header = TRUE)[,-1] 
sigSSTF3yr <- read.csv("output/uniSBI_sigSSTF3yr.csv",header = TRUE)[,-1] 

adjR2B_aggre <- matrix(0,nrow=28,ncol=5)

sig_Faggre <- matrix(0,nrow=28,ncol=5)
sig_SSTaggre <- matrix(0,nrow=28,ncol=5)
sig_SSTFaggre <- matrix(0,nrow=28,ncol=5)

for (i in 1:28) {
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

write.csv(adjR2B_aggre,"output/uniSBI_adjR2B_aggre.csv")
write.csv(sig_Faggre,"output/uniSBI_sigF_aggre.csv")
write.csv(sig_SSTaggre,"output/uniSBI_sigSST_aggre.csv")
write.csv(sig_SSTFaggre,"output/uniSBI_sigSSTF_aggre.csv")
