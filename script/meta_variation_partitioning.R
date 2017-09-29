# Meta-analysis of variation partitioning
# 2017/06/23 Crystal
# 2017/07/17 Merge fractionLM.R
# 2017/8/17 Use F/M instead of F in meanSST_F

setwd("~/Desktop/fishsize_varpart/script/")

lifehist_habitat <- read.csv("~/Desktop/fishsize_varpart/data/lifehist_habitat.csv",header = TRUE)
meanSST_F <- read.csv("~/Desktop/fishsize_varpart/output/meanSST_F.csv",header = TRUE)
cvSST_F <- read.csv("~/Desktop/fishsize_varpart/output/cvSST_F.csv",header = TRUE)

# combine the columns as variable
mean_cv <- cbind(meanSST_F[,-1],cvSST_F[,-1])
colnames(mean_cv)<-c("MeanF_M","MeanTemp","cvF_M","cvTemp")

# select the best based on total R-square
adjR2B <- read.csv("~/Desktop/fishsize_varpart/output/adjR2Bout.csv",header = TRUE)[,-1] 
adjR2B1yr <- read.csv("~/Desktop/fishsize_varpart/output/adjR2Bout1yr.csv",header = TRUE)[,-1]
adjR2B3yr <- read.csv("~/Desktop/fishsize_varpart/output/adjR2Bout3yr.csv",header = TRUE)[,-1]
colnames(adjR2B) <- c("Fishing","Interaction","Temperature","Total")
colnames(adjR2B1yr) <- c("Fishing","Interaction","Temperature","Total")
colnames(adjR2B3yr) <- c("Fishing","Interaction","Temperature","Total")

bestR2B <- max.col(cbind(adjR2B$Total,adjR2B1yr$Total,adjR2B3yr$Total))

sig <- read.csv("sigout.csv",header = TRUE)[,-1] 
sig1yr <- read.csv("sigout1yr.csv",header = TRUE)[,-1] 
sig3yr <- read.csv("sigout3yr.csv",header = TRUE)[,-1] 
colnames(sig) <- c("Fishing","Temperature")
colnames(sig1yr) <- c("Fishing","Temperature")
colnames(sig3yr) <- c("Fishing","Temperature")

adjR2B_aggre <- matrix(0,nrow=27,ncol=4)
sig_aggre <- matrix(0,nrow=27,ncol=2)

for (i in 1:27) {
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

output <- data.frame(lifehist_habitat,mean_cv,adjR2B,sig_aggre)

# Univariate linear regression on explained fraction of fishing/temperature v.s. Lifehist

coefSig <- matrix(NA,4,8)

#fishing
m_A50 <- lm(Fishing~1+A50,data=output)
coefSig[1,1] <- coef(m_A50)[2]
coefSig[2,1] <- anova(m_A50)$`Pr(>F)`[1]
m_L50 <- lm(Fishing~1+L50,data=output)
coefSig[1,2] <- coef(m_L50)[2]
coefSig[2,2] <- anova(m_L50)$`Pr(>F)`[1]
m_Linf <- lm(Fishing~1+Linf,data=output)
coefSig[1,3] <- coef(m_Linf)[2]
coefSig[2,3] <- anova(m_Linf)$`Pr(>F)`[1]
m_K <-lm(Fishing~1+K,data=output)
coefSig[1,4] <- coef(m_K)[2]
coefSig[2,4] <- anova(m_K)$`Pr(>F)`[1]
m_MeanF_M <- lm(Fishing~1+MeanF_M,data=output)
coefSig[1,5] <- coef(m_MeanF_M)[2]
coefSig[2,5] <- anova(m_MeanF_M)$`Pr(>F)`[1]
m_cvF_M <- lm(Fishing~1+cvF_M,data=output)
coefSig[1,6] <- coef(m_cvF_M)[2]
coefSig[2,6] <- anova(m_cvF_M)$`Pr(>F)`[1]
m_meanTemp <- lm(Fishing~1+MeanTemp,data=output)
coefSig[1,7] <- coef(m_meanTemp)[2]
coefSig[2,7] <- anova(m_meanTemp)$`Pr(>F)`[1]
m_cvTemp <- lm(Fishing~1+cvTemp,data=output)
coefSig[1,8] <- coef(m_cvTemp)[2]
coefSig[2,8] <- anova(m_cvTemp)$`Pr(>F)`[1]

#temperature
m_A50 <- lm(Temperature~1+A50,data=output)
coefSig[3,1] <- coef(m_A50)[2]
coefSig[4,1] <- anova(m_A50)$`Pr(>F)`[1]
m_L50 <- lm(Temperature~1+L50,data=output)
coefSig[3,2] <- coef(m_L50)[2]
coefSig[4,2] <- anova(m_L50)$`Pr(>F)`[1]
m_Linf <- lm(Temperature~1+Linf,data=output)
coefSig[3,3] <- coef(m_Linf)[2]
coefSig[4,3] <- anova(m_Linf)$`Pr(>F)`[1]
m_K <-lm(Temperature~1+K,data=output)
coefSig[3,4] <- coef(m_K)[2]
coefSig[4,4] <- anova(m_K)$`Pr(>F)`[1]
m_MeanF_M <- lm(Temperature~1+MeanF_M,data=output)
coefSig[3,5] <- coef(m_MeanF_M)[2]
coefSig[4,5] <- anova(m_MeanF_M)$`Pr(>F)`[1]
m_cvF_M <- lm(Temperature~1+cvF_M,data=output)
coefSig[3,6] <- coef(m_cvF_M)[2]
coefSig[4,6] <- anova(m_cvF_M)$`Pr(>F)`[1]
m_meanTemp <- lm(Temperature~1+MeanTemp,data=output)
coefSig[3,7] <- coef(m_meanTemp)[2]
coefSig[4,7] <- anova(m_meanTemp)$`Pr(>F)`[1]
m_cvTemp <- lm(Temperature~1+cvTemp,data=output)
coefSig[3,8] <- coef(m_cvTemp)[2]
coefSig[4,8] <- anova(m_cvTemp)$`Pr(>F)`[1]

coefSigOut <- data.frame(coefSig,row.names = c("coefF","sigF","coefT","sigT"))
colnames(coefSigOut) <- c("A50","L50","Linf","K","MeanF_M_M","cvF_M_M","meanTemp","cvTemp")

write.csv(output,"~/Desktop/fishsize_varpart/output/variationpartitioning_lifehist_mean_cv.csv")
write.csv(sig_aggre,"~/Desktop/fishsize_varpart/output/sigout_aggre.csv")
write.csv(coefSigOut,file = "~/Desktop/fishsize_varpart/output/fractionLMout.csv")
