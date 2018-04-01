# Univariate SBI for comparision
# Crystal Tu 2016/12/5
# 20170207 Add Atlantic
# 20170311 With revised North Sea length frequency 
# 20180214 With survey temperature and length composition from Alaska (in progress)
# 20180228 With time lag

# functions for indicators 
# (assume already run preset.r of all regions and run as part of do.r)
source("~/Desktop/fishsize_varpart/script/indicators.R")

# West US/Pacific 
sst <- read.table("~/Desktop/fishsize_varpart/data/WestUS/WestUSannualSSTfrom1948.txt",header = FALSE)
exploitation <- read.csv("~/Desktop/fishsize_varpart/data/FishingMortality.csv",header=TRUE)

lag=3

colnames(sst)<-c("YEAR","TEMP")
sigSST <- matrix(NA,28,5)
sigF <- matrix(NA,28,5)
sigSSTF <- matrix(NA,28,5)
r2B <- matrix(NA,28,5)  
year <- matrix(NA,28,2)

# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[1,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(12,80,2))
L95 <- L95perc(B,size)
mean <- meansize(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[1,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[1,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[1,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[1,j] <- summary(Breg)$adj.r.squared
}

# Chilipepper RF
B <- chilipepperRFfnum
idxF <- match(chilipepperRFf$YEAR,exploitation$YEAR)
idxSST <- match(chilipepperRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[2,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(16,52,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[2,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[2,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[2,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[2,j] <- summary(Breg)$adj.r.squared
}

#Darkblotched RF
B <- darkblotchedRFfnum
idxF <- match(darkblotchedRFf$YEAR,exploitation$YEAR)
idxSST <- match(darkblotchedRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[3,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,c(6:33,seq(35,51,2)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[3,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[3,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[3,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[3,j] <- summary(Breg)$adj.r.squared
}

#Dover sole
B <- dovernum
idxF <- match(doverS$YEAR,exploitation$YEAR)
idxSST <- match(doverS$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dover.Sole[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[4,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(22,60,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[4,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[4,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[4,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[4,j] <- summary(Breg)$adj.r.squared
}

#English sole
B <- englishFnum
idxF <- match(englishF$YEAR,exploitation$YEAR)
idxSST <- match(englishF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$English.Sole[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[5,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(11,45,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[5,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[5,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[5,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[5,j] <- summary(Breg)$adj.r.squared
}

#Lingcod
B <- lingcodFnum
idxF <- match(lingcodF$YEAR,exploitation$YEAR)
idxSST <- match(lingcodF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Lingcod[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[6,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(28,110,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[6,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[6,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[6,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[6,j] <- summary(Breg)$adj.r.squared
}

#Longspine
B <- longspineNum
idxF <- match(longspine$YEAR,exploitation$YEAR)
idxSST <- match(longspine$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[7,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(5,35,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[7,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[7,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[7,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[7,j] <- summary(Breg)$adj.r.squared
}

#Petrale Sole
B <- petraleSoleFnum
idxF <- match(petraleSole$YEAR,exploitation$YEAR)
idxSST <- match(petraleSole$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Petrale.sole[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[8,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(12,62,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[8,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[8,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[8,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[8,j] <- summary(Breg)$adj.r.squared
}

#Sardine
B <- sardineNum[-4,] #exclude 1984
idxF <- match(sardineS1$YEAR,exploitation$YEAR)
idxSST <- match(sardineS1$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Sardine[idxF],sst$TEMP[idxSST-lag]))[-4,]
colnames(A)<-c("F","SST")
year[9,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(9,26,0.5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[9,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[9,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[9,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[9,j] <- summary(Breg)$adj.r.squared
}

#Splitnose RF
B <- splitnoseRFfnum
idxF <- match(splitnoseRFf$YEAR,exploitation$YEAR)
idxSST <- match(splitnoseRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[10,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(2,56,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[10,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[10,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[10,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[10,j] <- summary(Breg)$adj.r.squared
}

#Yelloweye RF
B <- yelloweyeRFnum
idxF <- match(yelloweyeRF$YEAR,exploitation$YEAR)
idxSST <- match(yelloweyeRF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF],sst$TEMP[idxSST-lag]))
colnames(A)<-c("F","SST")
year[11,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(16,88,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[11,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[11,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[11,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[11,j] <- summary(Breg)$adj.r.squared
}

#Alaska
#AI pollock
newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[12,j] <- NA
  sigF[12,j] <- NA
  sigSSTF[12,j] <- NA
  r2B[12,j] <- NA
}

#GOA pollock
newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[13,j] <- NA
  sigF[13,j] <- NA
  sigSSTF[13,j] <- NA
  r2B[13,j] <- NA
}

#BSAI flathead sole
ebs_flatheadSoleAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/ebs_flatheadSoleAnnual.csv",header = TRUE)
ebsBottomT <- ebs_flatheadSoleAnnual[,c(1,5)]

B <- BSAIflatheadSoleNum
idxBottomT <- match(BSAIflatheadSoleF$YEAR,ebsBottomT$YEAR)
idxBottomT <- idxBottomT[is.finite(idxBottomT)]
idxF <- match(ebsBottomT$YEAR[idxBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]

idxF <- idxF[seq(1+lag,length(idxF))]
idxBottomT <- idxBottomT[seq(1,length(idxBottomT)-lag)]

A <- data.frame(cbind(exploitation$Flathead.sole.BSAI[idxF],ebsBottomT$BOT_TEMP[idxBottomT-1]))

idxB <- match(exploitation$YEAR[idxF],BSAIflatheadSoleF$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[14,] <- c(min(ebsBottomT$YEAR[idxBottomT]),max(ebsBottomT$YEAR[idxBottomT]))

size <- fishsize(B,c(seq(6,40,2),seq(43,58,3)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[14,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[14,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[14,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[14,j] <- summary(Breg)$adj.r.squared
}

#GOA flathead sole
newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[15,j] <- NA
  sigF[15,j] <- NA
  sigSSTF[15,j] <- NA
  r2B[15,j] <- NA
}

#EBS pacific cod
ebs_pacificCodAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/ebs_pacificCodAnnual.csv",header = TRUE)
ebsBottomT <- ebs_pacificCodAnnual[,c(1,5)]

B <- EBSpacifcCodNum

idxBottomT <- match(EBSpacifcCod$YEAR,ebsBottomT$YEAR)
idxBottomT <- idxBottomT[is.finite(idxBottomT)]
idxF <- match(ebsBottomT$YEAR[idxBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxF <- idxF[seq(1+lag,length(idxF))]
idxBottomT <- idxBottomT[seq(1,length(idxBottomT)-lag)]
A <- data.frame(cbind(exploitation$Pacific.cod.BSAI[idxF],ebsBottomT$BOT_TEMP[idxBottomT]))

idxB <- match(exploitation$YEAR[idxF],EBSpacifcCod$YEAR)
B <- B[idxB,]

colnames(A)<-c("F","SST")
year[16,] <- c(min(ebsBottomT$YEAR[idxBottomT]),max(sst$YEAR[idxBottomT]))

size <- fishsize(B,c(seq(9,42,3),seq(45,105,5)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[16,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[16,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[16,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[16,j] <- summary(Breg)$adj.r.squared
}


#GOA pacific cod
newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[17,j] <- NA
  sigF[17,j] <- NA
  sigSSTF[17,j] <- NA
  r2B[17,j] <- NA
}

#GOA rex sole
newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[18,j] <- NA
  sigF[18,j] <- NA
  sigSSTF[18,j] <- NA
  r2B[18,j] <- NA
}

#North Sea
YEAR<-1977:2014
#cod
idxF <- match(YEAR,cod$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(cod$Year[idxF],codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

year[19,]<-c(min(cod$Year[idxF]),max(cod$Year[idxF]))

size <- fishsize(B,seq(10,1400,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[19,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[19,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[19,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[19,j] <- summary(Breg)$adj.r.squared
}

#haddock
idxF <- match(YEAR,haddock$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(haddock$Year[idxF],haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

year[20,]<-c(min(haddock$Year[idxF]),max(haddock$Year[idxF]))

size <- fishsize(B,seq(10,870,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[20,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[20,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[20,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[20,j] <- summary(Breg)$adj.r.squared
}

#herring
idxF <- match(YEAR,herring$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(herring$Year[idxF],herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

year[21,]<-c(min(herring$Year[idxF]),max(herring$Year[idxF]))

size <- fishsize(B,seq(5,380,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[21,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[21,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[21,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[21,j] <- summary(Breg)$adj.r.squared
}

#mackerel
idxF <- match(YEAR,mackerel$Year)
idxF <- idxF[!is.na(idxF)]
idxbottomT <- match(mackerel$Year[idxF],NSavg_tempQ1$YEAR)

idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[!is.na(idxbottomT)]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(mackerel$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(mackerel$Year[idxF],mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

year[22,]<-c(min(mackerel$Year[idxF]),max(mackerel$Year[idxF]))

size <- fishsize(B,seq(10,580,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[22,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[22,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[22,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[22,j] <- summary(Breg)$adj.r.squared
}

#norwaypout
idxF <- match(YEAR,norwaypout$Year)
idxF <- idxF[!is.na(idxF)]
idxbottomT <- match(norwaypout$Year[idxF],NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[!is.na(idxbottomT)]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(norwaypout$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(norwaypout$Year[idxF[!is.na(idxF)]],norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

year[23,]<-c(min(norwaypout$Year[idxF]),max(norwaypout$Year[idxF]))

size <- fishsize(B,seq(10,320,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[23,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[23,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[23,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[23,j] <- summary(Breg)$adj.r.squared
}

#plaice
idxF <- match(YEAR,plaice$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(plaice$Year[idxF],plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

year[24,]<-c(min(plaice$Year[idxF]),max(plaice$Year[idxF]))

size <- fishsize(B,seq(10,650,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[24,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[24,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[24,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[24,j] <- summary(Breg)$adj.r.squared
}

#saithe
idxF <- match(YEAR,saithe$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(saithe$Year[idxF],saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

year[25,]<-c(min(saithe$Year[idxF]),max(saithe$Year[idxF]))

size <- fishsize(B,seq(10,1180,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[25,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[25,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[25,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[25,j] <- summary(Breg)$adj.r.squared
}

#sole
idxF <- match(YEAR,sole$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(sole$Year[idxF],soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

year[26,]<-c(min(sole$Year[idxF]),max(sole$Year[idxF]))

size <- fishsize(B,seq(10,510,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[26,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[26,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[26,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[26,j] <- summary(Breg)$adj.r.squared
}

#sprat
idxF <- match(YEAR,sprat$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")

idxB<- match(sprat$Year[idxF],spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

year[27,]<-c(min(sprat$Year[idxF]),max(sprat$Year[idxF]))

size <- fishsize(B,seq(5,245,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[27,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[27,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[27,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[27,j] <- summary(Breg)$adj.r.squared
}

#whiting
idxF <- match(YEAR,whiting$Year)
idxF <- idxF[!is.na(idxF)]
idxbottomT <- match(whiting$Year[idxF],NSavg_tempQ1$YEAR)

idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[!is.na(idxbottomT)]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(whiting$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","TEMP")
idxB<- match(whiting$Year[idxF],whitingLF$Year)

B <- whitingLF[idxB,-1] #LF w/o year

year[28,]<-c(min(whiting$Year[idxF]),max(whiting$Year[idxF]))

size <- fishsize(B,seq(10,690,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$TEMP+A$F+A$TEMP*A$F)
  sigSST[28,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[28,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[28,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[28,j] <- summary(Breg)$adj.r.squared
}

colnames(sigF) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSST) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSSTF) <- c("mean","L95","shannon","evenness","skewness")

write.csv(sigF,file = "output/uniSBI_sigF3yr.csv")
write.csv(sigSST,file = "output/uniSBI_sigSST3yr.csv")
write.csv(sigSSTF,file = "output/uniSBI_sigSSTF3yr.csv")
write.csv(r2B,file = "output/uniSBI_r2B3yr.csv")
write.csv(year,file="output/year_analyzed3yr.csv")