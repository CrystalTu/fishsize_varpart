# 20161210 Variation partitioning v.s. univariate SBIs
# 20170312 Plot with Atlantic species
# 20170712 Plot with new analysis (use cv instead of var, aggregate of 0, 1, 3 yr lag)
# 20180219 Plot with new analysis result (use survey data in Alaska)
# 20180228 Plot with time-lag result
library(ggplot2)
source("~/Desktop/fishsize_varpart/script/multiplot.R")
setwd("~/Desktop/fishsize_varpart/revision")

sigF<-read.csv("output/uniSBI_sigF_aggre.csv",header = TRUE)
sigSST<-read.csv("output/uniSBI_sigSST_aggre.csv",header = TRUE)
sigout<-read.csv("output/sigout_aggre.csv",header = TRUE)

# Exploitation
variationpartition <- c(rep(sigout$sigF,5))
univariateSBI <- c(sigF$mean,sigF$L95,sigF$shannon,sigF$evenness,sigF$skewness)
type <- factor(c(rep("mean",28),rep("L95",28),rep("shannon",28),rep("evenness",28),rep("skewness",28)))
exploitation <-data.frame(variationpartition,univariateSBI,type)

p1 <- ggplot(exploitation,aes(x=variationpartition,y=univariateSBI))
p1 <- p1+geom_point(aes(colour=factor(type)))+geom_abline(intercept = 0,slope = 1)
p1 <- p1+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_classic()
p1 <- p1+theme(panel.border = element_rect(colour = "black", fill=NA, size=1),legend.title=element_blank())
p1 <- p1+xlab("p values of variation partitioning")+ylab("p values of SBIs-based analysis")
p1 <- p1+ggtitle("(a) Fishing")+theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
#ggsave(file="comparisonExploitation.eps")
#ggsave(file="comparisonExploitation_revised.eps")

# Temperature
variationpartition <- c(rep(sigout$sigF,5))
univariateSBI <- c(sigSST$mean,sigSST$L95,sigSST$shannon,sigSST$evenness,sigSST$skewness)
type <- factor(c(rep("mean",28),rep("L95",28),rep("shannon",28),rep("evenness",28),rep("skewness",28)))
temperature <- data.frame(variationpartition,univariateSBI,type)

p2 <- ggplot(temperature,aes(x=variationpartition,y=univariateSBI))
p2 <- p2+geom_point(aes(colour=factor(type)))+geom_abline(intercept = 0,slope = 1)
p2 <- p2+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_classic()
p2 <- p2+theme(panel.border = element_rect(colour = "black", fill=NA, size=1),legend.title=element_blank())
p2 <- p2+xlab("p values of variation partitioning")+ylab("p values of SBIs-based analysis")
p2 <- p2+ggtitle("(b) Temperature")+theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))
#ggsave(file="comparisonTemperature.eps")
#ggsave(file="comparisonTemperature_revised.eps")

#binomial test
successF <- sum(exploitation$variationpartition < exploitation$univariateSBI,na.rm = TRUE)
successT <- sum(temperature$variationpartition < temperature$univariateSBI,na.rm = TRUE)
bi_successF <- c(binom.test(successF,28*5,0.5,"greater")$estimate,binom.test(successF,28*5,0.5,"greater")$p.value)
bi_successT <- c(binom.test(successT,28*5,0.5,"greater")$estimate,binom.test(successT,28*5,0.5,"greater")$p.value)

out <- matrix(c(1,2),nrow=1,byrow=TRUE)
setEPS()
postscript("output/Fig5_comparison.eps",width=13,height=6)
multiplot(p1,p2,cols=2)
dev.off()