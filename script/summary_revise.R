# Summary analysis & plot of variance partitioning
# 20161109 Crystal (1st ed.)
# 20170207 Add Atlantic
# 20170311 Plot Atlantic again with revised North Sea length frequency
# 20170421 Black&White plot
# 20170428 Black&White plot
# 20170517 Final plot
# 20170710 Plot with CV as index of variation (instead of variance in previous analysis)
# 20170719 Add stepwise forward selection for all variable
# 20170819 Use mortality ratio (F/M)
# 20180214 Revise with new result with survey temperature and length composition fro Alaska
# 20180226 Add time lag result
# 20180306 Corrected subplot annotation

library(lme4)
library(ggplot2)
library(Rcpp)
library(gridExtra)
library(dplyr)

source("~/Desktop/fishsize_varpart/script/multiplot.R")

result <- read.csv('output/variationPartitioning_lifehist_mean_cv_withlag.csv',header = TRUE)
nSp <- 28

# Rearrange the data in long format for overall boxplot
value <- c(result$Fishing,result$Interaction,result$Temperature)
habitat <- factor(rep(result$Habitat,3))
area <- factor(rep(result$Area2,3))
part <- factor(c(rep("Fishing",nSp),rep("Interaction",nSp),rep("Temperature",nSp)))
test <- data.frame(value,area,habitat,part)

b <- ggplot(test,aes(x=part,y=value))+geom_boxplot(show.legend = FALSE)
b <- b+theme_classic()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
b <- b+xlab("")+ylab("Fraction of variation explained")
b <- b+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
ggsave(file="output/Fig1_totalvariation.eps")

# pairwise comparison between fishing,interaction and temperature
mod1 <- anova(lm(value~part,data=test))
pair1 <- pairwise.t.test(test$value,test$part,p.adj = "bonf")

# Effect from area
p <- ggplot(test,aes(x = area, y = value, fill=part)) + theme(panel.background = element_rect(fill = "white", colour = "white"),axis.line=element_line(size=0.5,colour="black"))+ geom_boxplot() 
p <- p+theme_classic()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p <- p+xlab("")+ylab("Fraction of variation explained")
p <- p+theme(legend.title=element_blank(),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
ggsave(file="output/Fig2_area.eps",width=7,height=5.6,units = "in")

# pairwise comparison in each region
mod2 <- anova(lm(value~part*area,data=test))
pair2_1 <- pairwise.t.test(filter(test,area=="Alaska")$value,filter(test,area=="Alaska")$part,p.adj="bonf")
pair2_2 <- pairwise.t.test(filter(test,area=="North Sea")$value,filter(test,area=="North Sea")$part,p.adj="bonf")               
pair2_3 <- pairwise.t.test(filter(test,area=="West US")$value,filter(test,area=="West US")$part,p.adj="bonf")               

# Effect from habitat
theme_update()
p <- ggplot(test,aes(x = habitat, y = value, fill=part)) + theme(panel.background = element_rect(fill = "white", colour = "white"),axis.line=element_line(size=0.5,colour="black"))+ geom_boxplot() 
p <- p+theme_classic()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
p <- p+xlab("")+ylab("Fraction of variation explained")
p <- p+theme(legend.title=element_blank(),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
ggsave(file="output/Fig3_habitat.eps",width=7,height=5.6,units = "in")

mod3 <- anova(lm(value~part*habitat,data=test))
pair3_1 <- pairwise.t.test(filter(test,habitat=="bathydemersal")$value,filter(test,habitat=="bathydemersal")$part,p.adj="bonf")
pair3_2 <- pairwise.t.test(filter(test,habitat=="benthopelagic")$value,filter(test,habitat=="benthopelagic")$part,p.adj="bonf")
pair3_3 <- pairwise.t.test(filter(test,habitat=="demersal")$value,filter(test,habitat=="demersal")$part,p.adj="bonf")
pair3_4 <- pairwise.t.test(filter(test,habitat=="pelagic")$value,filter(test,habitat=="pelagic")$part,p.adj="bonf")

# LMM on explained variance v.s. life history traits
library(lmerTest)
# Fishing
m <- lmer(Fishing~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_K <- lmer(Fishing~K+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_K)))
m.p= anova(m_K,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
b=coefs[1,1]
a=coefs[2,1]
p1 <- ggplot(result,aes(x=K, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5) 
p1 <- p1 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$K),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p1 <- p1 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p1 <- p1 + ggtitle("(a)") + ylab("Fraction of variation explained by fishing") 
#ggsave(file="GLMM_K.eps")

m_Linf <- lmer(Fishing~Linf+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_Linf)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_Linf,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p2 <- ggplot(result,aes(x=Linf, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5)
p2 <- p2 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$Linf),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p2 <- p2 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p2 <- p2 + ggtitle("(b)") + ylab("Fraction of variation explained by fishing") 
#ggsave(file="GLMM_Linf.eps")

m_A50 <- lmer(Fishing~A50+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_A50)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_A50,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p3 <- ggplot(result,aes(x=A50, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5)
p3 <- p3 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$A50),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p3 <- p3 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p3 <- p3 + ggtitle("(c)") + ylab("Fraction of variation explained by fishing") 
#ggsave(file="GLMM_A50.eps")

m_L50 <- lmer(Fishing~L50+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_L50)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_L50,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p4 <- ggplot(result,aes(x=L50, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5)
p4 <- p4 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$L50),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p4 <- p4 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p4 <- p4 + ggtitle("(d)") + ylab("Fraction of variation explained by fishing") 
#ggsave(file="GLMM_L50.eps")

out <- matrix(c(1,2,3,4),nrow=2,byrow=TRUE)
setEPS()
postscript("output/FigS4_LMM_Fishing_lifehistory_final.eps",width=15,height=12)
multiplot(p1,p2,p3,p4,layout=out)
dev.off()

#Fishing v.s. mean/cv of temperature and mortality ratio
m <- lmer(Fishing~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_MeanTemp <- lmer(Fishing~MeanTemp+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_MeanTemp)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_MeanTemp,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p5 <- ggplot(result,aes(x=MeanTemp, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5)
p5 <- p5 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$MeanTemp),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p5 <- p5 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p5 <- p5 + ggtitle("(a)")+ ylab("Fraction of variation explained by fishing") 
#ggsave(file="GLMM_MeanTemp.eps")

m <- lmer(Fishing~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_cvTemp <- lmer(Fishing~cvTemp+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_cvTemp)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_cvTemp,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p6 <- ggplot(result,aes(x=cvTemp, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5)
p6 <- p6 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$cvTemp),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p6 <- p6 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p6 <- p6 + ggtitle("(b)")+ ylab("Fraction of variation explained by fishing") 
##ggsave(file="GLMM_cvTemp.eps")

m <- lmer(Fishing~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_MeanF_M <- lmer(Fishing~MeanF_M+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_MeanF_M)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_MeanF_M,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p7 <- ggplot(result,aes(x=MeanF_M, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5)
p7 <- p7 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$MeanF_M),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p7 <- p7 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p7 <- p7 + ggtitle("(c)")+ ylab("Fraction of variation explained by fishing") 
#ggsave(file="GLMM_MeanF.eps")

m <- lmer(Fishing~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_cvF_M <- lmer(Fishing~cvF_M+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_cvF_M)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_cvF_M,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p8 <- ggplot(result,aes(x=cvF_M, y=Fishing)) + geom_point(aes(shape=factor(Habitat)),size=5)
p8 <- p8 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$cvF_M),y=max(result$Fishing),hjust=1,vjust=1,size=8,label=p_text)
p8 <- p8 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p8 <- p8 + ggtitle("(d)")+ ylab("Fraction of variation explained by fishing") 
#ggsave(file="GLMM_cvF.eps")

out <- matrix(c(1,2,3,4),nrow=2,byrow=TRUE)
setEPS()
postscript("output/FigS5_LMM_Fishing_FandT_naturalmortality_final.eps",width=15,height=12)
multiplot(p5,p6,p7,p8,layout=out)
dev.off()

# Temperature
m <- lmer(Temperature~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_K <- lmer(Temperature~K+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_K)))
m.p= anova(m_K,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
b=coefs[1,1]
a=coefs[2,1]
p1 <- ggplot(result,aes(x=K, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5)
p1 <- p1 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$K),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p1 <- p1 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p1 <- p1 + ggtitle("(a)")+ ylab("Fraction of variation explained by temperature") 
#ggsave(file="GLMM_TK.eps")

m_Linf <- lmer(Temperature~Linf+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_Linf)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_Linf,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p2 <- ggplot(result,aes(x=Linf, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5) 
p2 <- p2 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$Linf),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p2 <- p2 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p2 <- p2 + ggtitle("(b)")+ ylab("Fraction of variation explained by temperature") 
#ggsave(file="GLMM_TLinf.eps")

m_A50 <- lmer(Temperature~A50+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_A50)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_A50,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p3 <- ggplot(result,aes(x=A50, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5)
p3 <- p3 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$A50),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p3 <- p3 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p3 <- p3 + ggtitle("(c)")+ ylab("Fraction of variation explained by temperature") 
#ggsave(file="GLMM_TA50.eps") 

m <- lmer(Temperature~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_L50 <- lmer(Temperature~L50+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_L50)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_L50,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p4 <- ggplot(result,aes(x=L50, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5) 
p4 <- p4 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$L50),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p4 <- p4 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p4 <- p4 + ggtitle("(d)")+ ylab("Fraction of variation explained by temperature") 
#ggsave(file="GLMM_TL50.eps")

out <- matrix(c(1,2,3,4),nrow=2,byrow=TRUE)
setEPS()
postscript("output/FigS6_LMM_Temperature_lifehistory_final.eps",width=15,height=12)
multiplot(p1,p2,p3,p4,layout=out)
dev.off()

# Temperature v.s mean/cv of temeprature and mortality ratio
m <- lmer(Temperature~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_MeanTemp <- lmer(Temperature~MeanTemp+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_MeanTemp)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_MeanTemp,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p5 <- ggplot(result,aes(x=MeanTemp, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5) 
p5 <- p5 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$MeanTemp),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p5 <- p5 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p5 <- p5 + ggtitle("(a)")+ ylab("Fraction of variation explained by temperature") 

m_cvTemp <- lmer(Temperature~cvTemp+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_cvTemp)))
m.p= anova(m_cvTemp,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
b=coefs[1,1]
a=coefs[2,1]
p6 <- ggplot(result,aes(x=cvTemp, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5) 
p6 <- p6 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$cvTemp),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p6 <- p6 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p6 <- p6 + ggtitle("(b)")+ ylab("Fraction of variation explained by temperature") 

m <- lmer(Temperature~1+(1|Habitat),data=result,REML=FALSE) #reduced model
m_MeanF_M <- lmer(Temperature~MeanF_M+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_MeanF_M)))
b=coefs[1,1]
a=coefs[2,1]
m.p= anova(m_MeanF_M,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
p7 <- ggplot(result,aes(x=MeanF_M, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5)
p7 <- p7 + geom_abline(intercept = b, slope = a,linetype = "dashed") + annotate("text",x=max(result$MeanF_M),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p7 <- p7 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p7 <- p7 + ggtitle("(c)")+ ylab("Fraction of variation explained by temperature") 

m_cvF_M <- lmer(Temperature~cvF_M+(1|Habitat),data=result,REML=FALSE)
coefs <- data.frame(coef(summary(m_cvF_M)))
m.p= anova(m_cvF_M,m)
p_text = paste("p=",as.character(round(m.p$`Pr(>Chisq)`[2],digits = 3)))
b=coefs[1,1]
a=coefs[2,1]
p8 <- ggplot(result,aes(x=cvF_M, y=Temperature)) + geom_point(aes(shape=factor(Habitat)),size=5) 
p8 <- p8 + geom_abline(intercept = b, slope = a,linetype = "solid") + annotate("text",x=max(result$cvF_M),y=max(result$Temperature),hjust=1,vjust=1,size=8,label=p_text)
p8 <- p8 + theme_classic() + theme(text = element_text(size=16),panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p8 <- p8 + ggtitle("(d)")+ ylab("Fraction of variation explained by temperature") 

out <- matrix(c(1,2,3,4),nrow=2,byrow=TRUE)
setEPS()
postscript("output/Fig4_LMM_Temperature_FandT_naturalmortality_final.eps",width=15,height=12)
multiplot(p5,p6,p7,p8,layout=out)
dev.off()

# Univariate linear regression on explained fraction of fishing/temperature v.s. Lifehist
coefSig <- matrix(NA,4,8)

#fishing
m_A50 <- lm(Fishing~1+A50,data=result)
coefSig[1,1] <- coef(m_A50)[2]
coefSig[2,1] <- anova(m_A50)$`Pr(>F)`[1]
m_L50 <- lm(Fishing~1+L50,data=result)
coefSig[1,2] <- coef(m_L50)[2]
coefSig[2,2] <- anova(m_L50)$`Pr(>F)`[1]
m_Linf <- lm(Fishing~1+Linf,data=result)
coefSig[1,3] <- coef(m_Linf)[2]
coefSig[2,3] <- anova(m_Linf)$`Pr(>F)`[1]
m_K <-lm(Fishing~1+K,data=result)
coefSig[1,4] <- coef(m_K)[2]
coefSig[2,4] <- anova(m_K)$`Pr(>F)`[1]
m_MeanF_M <- lm(Fishing~1+MeanF_M,data=result)
coefSig[1,5] <- coef(m_MeanF_M)[2]
coefSig[2,5] <- anova(m_MeanF_M)$`Pr(>F)`[1]
m_cvF_M <- lm(Fishing~1+cvF_M,data=result)
coefSig[1,6] <- coef(m_cvF_M)[2]
coefSig[2,6] <- anova(m_cvF_M)$`Pr(>F)`[1]
m_meanTemp <- lm(Fishing~1+MeanTemp,data=result)
coefSig[1,7] <- coef(m_meanTemp)[2]
coefSig[2,7] <- anova(m_meanTemp)$`Pr(>F)`[1]
m_cvTemp <- lm(Fishing~1+cvTemp,data=result)
coefSig[1,8] <- coef(m_cvTemp)[2]
coefSig[2,8] <- anova(m_cvTemp)$`Pr(>F)`[1]

#temperature
m_A50 <- lm(Temperature~1+A50,data=result)
coefSig[3,1] <- coef(m_A50)[2]
coefSig[4,1] <- anova(m_A50)$`Pr(>F)`[1]
m_L50 <- lm(Temperature~1+L50,data=result)
coefSig[3,2] <- coef(m_L50)[2]
coefSig[4,2] <- anova(m_L50)$`Pr(>F)`[1]
m_Linf <- lm(Temperature~1+Linf,data=result)
coefSig[3,3] <- coef(m_Linf)[2]
coefSig[4,3] <- anova(m_Linf)$`Pr(>F)`[1]
m_K <-lm(Temperature~1+K,data=result)
coefSig[3,4] <- coef(m_K)[2]
coefSig[4,4] <- anova(m_K)$`Pr(>F)`[1]
m_MeanF_M <- lm(Temperature~1+MeanF_M,data=result)
coefSig[3,5] <- coef(m_MeanF_M)[2]
coefSig[4,5] <- anova(m_MeanF_M)$`Pr(>F)`[1]
m_cvF_M <- lm(Temperature~1+cvF_M,data=result)
coefSig[3,6] <- coef(m_cvF_M)[2]
coefSig[4,6] <- anova(m_cvF_M)$`Pr(>F)`[1]
m_meanTemp <- lm(Temperature~1+MeanTemp,data=result)
coefSig[3,7] <- coef(m_meanTemp)[2]
coefSig[4,7] <- anova(m_meanTemp)$`Pr(>F)`[1]
m_cvTemp <- lm(Temperature~1+cvTemp,data=result)
coefSig[3,8] <- coef(m_cvTemp)[2]
coefSig[4,8] <- anova(m_cvTemp)$`Pr(>F)`[1]

coefSigOut <- data.frame(coefSig,row.names = c("coefF","sigF","coefT","sigT"))
colnames(coefSigOut) <- c("A50","L50","Linf","K","MeanF_M","cvF_M","meanTemp","cvTemp")
write.csv(coefSigOut,file = "output/Table3_fractionLMout.csv")