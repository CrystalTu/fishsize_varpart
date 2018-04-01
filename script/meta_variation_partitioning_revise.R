# Meta-analysis of variation partitioning
# 2017/06/23 Crystal
# 2017/07/17 Merge fractionLM.R
# 2017/08/17 Use F/M instead of F in meanSST_F
# 2018/02/14 New result w/ survey bottom temperature and length composition
# 2018/02/28 Add time lag

lifehist_habitat <- read.csv("~/Desktop/fishsize_varpart/data/lifehist_habitat.csv",header = TRUE)
meanSST_F <- read.csv("output/meanSST_F.csv",header = TRUE)
cvSST_F <- read.csv("output/cvSST_F.csv",header = TRUE)

# combine the columns as variable
mean_cv <- cbind(meanSST_F[,-1],cvSST_F[,-1])
colnames(mean_cv)<-c("MeanF_M","MeanTemp","cvF_M","cvTemp")

# select the best based on total R-square
adjR2B <- read.csv("output/adjR2Bout.csv",header = TRUE)[,-1] 
adjR2B1yr <- read.csv("output/adjR2Bout1yrTemp.csv",header = TRUE)[,-1]
adjR2B3yr <- read.csv("output/adjR2Bout3yrTemp.csv",header = TRUE)[,-1]
adjR2B3yr[is.na(adjR2B3yr)] <- 0 #replace NA with 0
  
colnames(adjR2B) <- c("Fishing","Interaction","Temperature","Total")
colnames(adjR2B1yr) <- c("Fishing","Interaction","Temperature","Total")
colnames(adjR2B3yr) <- c("Fishing","Interaction","Temperature","Total")
bestR2B <- max.col(cbind(adjR2B$Total,adjR2B1yr$Total,adjR2B3yr$Total), ties.method = "first")
  
sig <- read.csv("output/sigout.csv",header = TRUE)[,-1] 
  sig1yr <- read.csv("output/sigout1yrTemp.csv",header = TRUE)[,-1] 
  sig3yr <- read.csv("output/sigout3yrTemp.csv",header = TRUE)[,-1] 
  colnames(sig) <- c("Fishing","Temperature")
  colnames(sig1yr) <- c("Fishing","Temperature")
  colnames(sig3yr) <- c("Fishing","Temperature")

adjR2B_aggre <- matrix(0,nrow=28,ncol=4)
sig_aggre <- matrix(0,nrow=28,ncol=2)

for (i in 1:28) {
  if (bestR2B[i] == 1) {
  adjR2B_aggre[i,] <- as.matrix(adjR2B[i,])[1:4]
  sig_aggre[i,] <- as.matrix(sig[i,])[1:2] }
  if (bestR2B[i] == 2) {
  adjR2B_aggre[i,] <- as.matrix(adjR2B1yr[i,])[1:4]
  sig_aggre[i,] <- as.matrix(sig1yr[i,])[1:2] }
  if (bestR2B[i] == 3) {
  adjR2B_aggre[i,] <- as.matrix(adjR2B3yr[i,])[1:4]
  sig_aggre[i,] <- as.matrix(sig3yr[i,])[1:2] }
}

colnames(adjR2B_aggre) <- c("Fishing","Interaction","Temperature","Total")
colnames(sig_aggre) <- c("sigF","sigT")
colnames(mean_cv)<-c("MeanF_M","MeanTemp","cvF_M","cvTemp")

output <- data.frame(lifehist_habitat,mean_cv,adjR2B_aggre,sig_aggre)
write.csv(sig_aggre,"output/sigout_aggre.csv",row.names = FALSE)
write.csv(output,"output/variationpartitioning_lifehist_mean_cv_withlag.csv")
write.csv(bestR2B,"output/bestR2B.csv")
