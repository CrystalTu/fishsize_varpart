# Univariate SBI for comparision
# Crystal Tu 2016/12/5
# 20170207 Add Atlantic
# 20170311 With revised North Sea length frequency 

# run preset.R in length/ first
source("~/Desktop/fishsize_varpart/script/indicators.R")
setwd("~/Desktop/fishsize_varpart/")

# West US/Pacific 
source("~/Desktop/fishsize_varpart/data/length/preset.R") #length 
sst <- read.table("~/Desktop/fishsize_varpart/data/WestUS/WestUSannualSSTfrom1948.txt",header = FALSE)
exploitation <- read.csv("~/Desktop/fishsize_varpart/data/WestUS/exploitation.csv",header=TRUE)

colnames(sst)<-c("YEAR","SST")
sigSST <- matrix(NA,27,5)
sigF <- matrix(NA,27,5)
sigSSTF <- matrix(NA,27,5)
r2B <- matrix(NA,27,5)  
year <- matrix(NA,27,2)

# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Dover.Sole[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$English.Sole[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Lingcod[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Petrale.sole[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Sardine[idxF],sst$SST[idxSST]))[-4,]
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
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF],sst$SST[idxSST]))
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
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF],sst$SST[idxSST]))
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
source("~/Desktop/fishsize_varpart/data/Alaska/preset.R")
setwd("~/Desktop/fishsize_varpart/")

#AI pollock
B <- AIpollocknum
idxF <- match(AIpollockLF$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxSST <- match(exploitation$YEAR[idxF],AIsst$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..AI.[idxF],AIsst$SST[idxSST]))
idxB <- match(exploitation$YEAR[idxF],AIpollockLF$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[12,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(10,104,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[12,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[12,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[12,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[12,j] <- summary(Breg)$adj.r.squared
}

#GOA pollock
B <- GOApollockNum
idxF <- match(GOApollock$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..GOA.[idxF],GAK1bottomT$T[idxGAK1]))
idxB <- match(exploitation$YEAR[idxF],GOApollock$YEAR)
B <- B[idxB,]
year[13,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

colnames(A)<-c("F","SST")
size <- fishsize(B,seq(5,75,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[13,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[13,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[13,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[13,j] <- summary(Breg)$adj.r.squared
}

#BSAI flathead sole
B <- BSAIflatheadSoleNum
idxsummerBottomT <- match(BSAIflatheadSoleF$YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Flathead.sole..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT-1]))

idxB <- match(exploitation$YEAR[idxF],BSAIflatheadSoleF$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[14,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

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
B <- GOAflatheadSoleNum
idxF <- match(GOAflatheadSole$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Flathead.sole..GOA.[idxF],GAK1bottomT$T[idxGAK1]))

idxB <- match(exploitation$YEAR[idxF],GOAflatheadSole$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[15,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(14,48,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[15,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[15,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[15,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[15,j] <- summary(Breg)$adj.r.squared
}

#BSAI pacific cod
B <- BSAIpacifcCodNum
idxsummerBottomT <- match(BSAIpacifcCod$YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Pacific.cod..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

idxB <- match(exploitation$YEAR[idxF],BSAIpacifcCod$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[16,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

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

#GOA rex sole
B <- GOArexSoleNum
idxF <- match(GOArexSoleF$YEAR,exploitation$YEAR)
idxGAK1 <- match(GOArexSoleF$YEAR,GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Rex.sole..GOA.[idxF],GAK1bottomT$T[idxGAK1]))
colnames(A)<-c("F","SST")
year[17,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(9,41,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST[17,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[17,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[17,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[17,j] <- summary(Breg)$adj.r.squared
}

#Atlantic
setwd("~/Desktop/fishsize_varpart/")
source("~/Desktop/fishsize_varpart/data/NorthSea/IBTSsurvey/preset.R") # get length frequency data
source("~/Desktop/fishsize_varpart/data/NorthSea/stockassessment/stockassessment.R") # get stockassessment result table
source("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/temperature.R") # get bottom T
setwd("~/Desktop/fishsize_varpart/")

YEAR<-1977:2014

#cod
idxF <- match(YEAR,cod$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

year[18,]<-c(min(cod$Year[idxF]),max(cod$Year[idxF]))

size <- fishsize(B,seq(10,1400,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[18,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[18,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[18,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[18,j] <- summary(Breg)$adj.r.squared
}

#haddock
idxF <- match(YEAR,haddock$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

year[19,]<-c(min(haddock$Year[idxF]),max(haddock$Year[idxF]))

size <- fishsize(B,seq(10,870,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[19,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[19,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[19,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[19,j] <- summary(Breg)$adj.r.squared
}

#herring
idxF <- match(YEAR,herring$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

year[20,]<-c(min(herring$Year[idxF]),max(herring$Year[idxF]))

size <- fishsize(B,seq(5,380,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[20,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[20,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[20,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[20,j] <- summary(Breg)$adj.r.squared
}

#mackerel
idxF <- match(YEAR,mackerel$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(mackerel$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

year[21,]<-c(min(mackerel$Year[idxF]),max(mackerel$Year[idxF]))

size <- fishsize(B,seq(10,580,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[21,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[21,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[21,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[21,j] <- summary(Breg)$adj.r.squared
}

#norwaypout
idxF <- match(YEAR,norwaypout$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(norwaypout$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

year[22,]<-c(min(norwaypout$Year[idxF]),max(norwaypout$Year[idxF]))

size <- fishsize(B,seq(10,320,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[22,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[22,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[22,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[22,j] <- summary(Breg)$adj.r.squared
}

#plaice
idxF <- match(YEAR,plaice$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

year[23,]<-c(min(plaice$Year[idxF]),max(plaice$Year[idxF]))

size <- fishsize(B,seq(10,650,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[23,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[23,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[23,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[23,j] <- summary(Breg)$adj.r.squared
}

#saithe
idxF <- match(YEAR,saithe$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

year[24,]<-c(min(saithe$Year[idxF]),max(saithe$Year[idxF]))

size <- fishsize(B,seq(10,1180,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[24,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[24,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[24,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[24,j] <- summary(Breg)$adj.r.squared
}

#sole
idxF <- match(YEAR,sole$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

year[25,]<-c(min(sole$Year[idxF]),max(sole$Year[idxF]))

size <- fishsize(B,seq(10,510,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[25,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[25,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[25,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[25,j] <- summary(Breg)$adj.r.squared
}

#sprat
idxF <- match(YEAR,sprat$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

year[26,]<-c(min(sprat$Year[idxF]),max(sprat$Year[idxF]))

size <- fishsize(B,seq(5,245,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[26,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[26,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[26,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[26,j] <- summary(Breg)$adj.r.squared
}

#whiting
idxF <- match(YEAR,whiting$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,spratLF$Year)
B <- whitingLF[idxB,-1] #LF w/o year

year[27,]<-c(min(whiting$Year[idxF]),max(whiting$Year[idxF]))

size <- fishsize(B,seq(10,690,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST[27,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF[27,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF[27,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B[27,j] <- summary(Breg)$adj.r.squared
}

colnames(sigF) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSST) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSSTF) <- c("mean","L95","shannon","evenness","skewness")

write.csv(sigF,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigF.csv")
write.csv(sigSST,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigSST.csv")
write.csv(sigSSTF,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigSSTF.csv")
write.csv(r2B,file = "~/Desktop/fishsize_varpart/output/uniSBI_r2B.csv")
write.csv(year,file="~/Desktop/fishsize_varpart/output/year_analyzed.csv")

# 1-yr lag
sigSST_1yr <- matrix(NA,27,5)
sigF_1yr <- matrix(NA,27,5)
sigSSTF_1yr <- matrix(NA,27,5)
r2B_1yr <- matrix(NA,27,5)  
year <- matrix(NA,27,2)

# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[1,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(12,80,2))
L95 <- L95perc(B,size)
mean <- meansize(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[1,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[1,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[1,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[1,j] <- summary(Breg)$adj.r.squared
}

# Chilipepper RF
B <- chilipepperRFfnum
idxF <- match(chilipepperRFf$YEAR,exploitation$YEAR)
idxSST <- match(chilipepperRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[2,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(16,52,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[2,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[2,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[2,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[2,j] <- summary(Breg)$adj.r.squared
}

#Darkblotched RF
B <- darkblotchedRFfnum
idxF <- match(darkblotchedRFf$YEAR,exploitation$YEAR)
idxSST <- match(darkblotchedRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[3,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,c(6:33,seq(35,51,2)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[3,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[3,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[3,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[3,j] <- summary(Breg)$adj.r.squared
}

#Dover sole
B <- dovernum
idxF <- match(doverS$YEAR,exploitation$YEAR)
idxSST <- match(doverS$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dover.Sole[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[4,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(22,60,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[4,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[4,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[4,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[4,j] <- summary(Breg)$adj.r.squared
}

#English sole
B <- englishFnum
idxF <- match(englishF$YEAR,exploitation$YEAR)
idxSST <- match(englishF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$English.Sole[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[5,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(11,45,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[5,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[5,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[5,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[5,j] <- summary(Breg)$adj.r.squared
}

#Lingcod
B <- lingcodFnum
idxF <- match(lingcodF$YEAR,exploitation$YEAR)
idxSST <- match(lingcodF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Lingcod[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[6,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(28,110,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[6,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[6,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[6,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[6,j] <- summary(Breg)$adj.r.squared
}

#Longspine
B <- longspineNum
idxF <- match(longspine$YEAR,exploitation$YEAR)
idxSST <- match(longspine$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[7,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(5,35,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[7,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[7,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[7,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[7,j] <- summary(Breg)$adj.r.squared
}

#Petrale Sole
B <- petraleSoleFnum
idxF <- match(petraleSole$YEAR,exploitation$YEAR)
idxSST <- match(petraleSole$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Petrale.sole[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[8,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(12,62,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[8,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[8,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[8,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[8,j] <- summary(Breg)$adj.r.squared
}

#Sardine
B <- sardineNum[-4,] #exclude 1984
idxF <- match(sardineS1$YEAR,exploitation$YEAR)
idxSST <- match(sardineS1$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Sardine[idxF-1],sst$SST[idxSST-1]))[-4,]
colnames(A)<-c("F","SST")
year[9,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(9,26,0.5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[9,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[9,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[9,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[9,j] <- summary(Breg)$adj.r.squared
}

#Splitnose RF
B <- splitnoseRFfnum
idxF <- match(splitnoseRFf$YEAR,exploitation$YEAR)
idxSST <- match(splitnoseRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[10,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(2,56,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[10,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[10,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[10,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[10,j] <- summary(Breg)$adj.r.squared
}

#Yelloweye RF
B <- yelloweyeRFnum
idxF <- match(yelloweyeRF$YEAR,exploitation$YEAR)
idxSST <- match(yelloweyeRF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF-1],sst$SST[idxSST-1]))
colnames(A)<-c("F","SST")
year[11,] <- c(min(sst$YEAR[idxSST-1]),max(sst$YEAR[idxSST-1]))

size <- fishsize(B,seq(16,88,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[11,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[11,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[11,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[11,j] <- summary(Breg)$adj.r.squared
}

#AI pollock
B <- AIpollocknum
idxF <- match(AIpollockLF$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxSST <- match(exploitation$YEAR[idxF],AIsst$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..AI.[idxF-1],AIsst$SST[idxSST-1]))
idxB <- match(exploitation$YEAR[idxF],AIpollockLF$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[12,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(10,104,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[12,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[12,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[12,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[12,j] <- summary(Breg)$adj.r.squared
}

#GOA pollock
B <- GOApollockNum
idxF <- match(GOApollock$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..GOA.[idxF-1],GAK1bottomT$T[idxGAK1-1]))
idxB <- match(exploitation$YEAR[idxF],GOApollock$YEAR)
B <- B[idxB,]
year[13,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

colnames(A)<-c("F","SST")
size <- fishsize(B,seq(5,75,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[13,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[13,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[13,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[13,j] <- summary(Breg)$adj.r.squared
}

#BSAI flathead sole
B <- BSAIflatheadSoleNum
idxsummerBottomT <- match(BSAIflatheadSoleF$YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Flathead.sole..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT-1]))

idxB <- match(exploitation$YEAR[idxF],BSAIflatheadSoleF$YEAR)
B <- B[idxB+1,]  #shift the length distribution instead of temperature and F due to data availability
colnames(A)<-c("F","SST")
year[14,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,c(seq(6,40,2),seq(43,58,3)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[14,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[14,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[14,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[14,j] <- summary(Breg)$adj.r.squared
}

#GOA flathead sole (DONE)
B <- GOAflatheadSoleNum
idxF <- match(GOAflatheadSole$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Flathead.sole..GOA.[idxF-1],GAK1bottomT$T[idxGAK1-1]))

idxB <- match(exploitation$YEAR[idxF],GOAflatheadSole$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[15,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(14,48,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[15,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[15,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[15,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[15,j] <- summary(Breg)$adj.r.squared
}

#BSAI pacific cod
YEAR<-seq(min(summerBottomT$YEAR),max(BSAIpacifcCod$YEAR)-1)
B <- BSAIpacifcCodNum
idxsummerBottomT <- match(YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Pacific.cod..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

YEAR<-seq(min(summerBottomT$YEAR)+1,max(BSAIpacifcCod$YEAR))
idxB <- match(YEAR,BSAIpacifcCod$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[16,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,c(seq(9,42,3),seq(45,105,5)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[16,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[16,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[16,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[16,j] <- summary(Breg)$adj.r.squared
}

#GOA rex sole
B <- GOArexSoleNum
idxF <- match(GOArexSoleF$YEAR,exploitation$YEAR)
idxGAK1 <- match(GOArexSoleF$YEAR,GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Rex.sole..GOA.[idxF-1],GAK1bottomT$T[idxGAK1-1]))
colnames(A)<-c("F","SST")
year[17,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(9,41,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_1yr[17,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[17,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[17,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[17,j] <- summary(Breg)$adj.r.squared
}

#Atlantic
YEAR<-1977:2014

#cod
idxF <- match(YEAR,cod$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,codLF$Year)
B <- codLF[idxB+1,-1] #LF w/o year

year[18,]<-c(min(cod$Year[idxF]),max(cod$Year[idxF]))

size <- fishsize(B,seq(10,1400,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[18,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[18,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[18,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[18,j] <- summary(Breg)$adj.r.squared
}

#haddock
idxF <- match(YEAR,haddock$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,haddockLF$Year)
B <- haddockLF[idxB+1,-1] #LF w/o year

year[19,]<-c(min(haddock$Year[idxF]),max(haddock$Year[idxF]))

size <- fishsize(B,seq(10,870,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[19,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[19,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[19,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[19,j] <- summary(Breg)$adj.r.squared
}

#herring
idxF <- match(YEAR,herring$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,herringLF$Year)
B <- herringLF[idxB+1,-1] #LF w/o year

year[20,]<-c(min(herring$Year[idxF-1]),max(herring$Year[idxF-1]))

size <- fishsize(B,seq(5,380,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[20,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[20,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[20,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[20,j] <- summary(Breg)$adj.r.squared
}

#mackerel
idxF <- match(YEAR,mackerel$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(mackerel$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,mackerelLF$Year)
B <- mackerelLF[idxB+1,-1] #LF w/o year

year[21,]<-c(min(mackerel$Year[idxF]),max(mackerel$Year[idxF]))

size <- fishsize(B,seq(10,580,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[21,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[21,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[21,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[21,j] <- summary(Breg)$adj.r.squared
}

#norwaypout
idxF <- match(YEAR,norwaypout$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(norwaypout$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,norwaypoutLF$Year)
B <- norwaypoutLF[idxB+1,-1] #LF w/o year

year[22,]<-c(min(norwaypout$Year[idxF]),max(norwaypout$Year[idxF]))

size <- fishsize(B,seq(10,320,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[22,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[22,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[22,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[22,j] <- summary(Breg)$adj.r.squared
}

#plaice
idxF <- match(YEAR,plaice$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,plaiceLF$Year)
B <- plaiceLF[idxB+1,-1] #LF w/o year

year[23,]<-c(min(plaice$Year[idxF]),max(plaice$Year[idxF]))

size <- fishsize(B,seq(10,650,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[23,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[23,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[23,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[23,j] <- summary(Breg)$adj.r.squared
}

#saithe
idxF <- match(YEAR,saithe$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,saitheLF$Year)
B <- saitheLF[idxB+1,-1] #LF w/o year

year[24,]<-c(min(saithe$Year[idxF]),max(saithe$Year[idxF]))

size <- fishsize(B,seq(10,1180,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[24,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[24,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[24,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[24,j] <- summary(Breg)$adj.r.squared
}

#sole
idxF <- match(YEAR,sole$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,soleLF$Year)
B <- soleLF[idxB+1,-1] #LF w/o year

year[25,]<-c(min(sole$Year[idxF]),max(sole$Year[idxF]))

size <- fishsize(B,seq(10,510,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[25,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[25,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[25,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[25,j] <- summary(Breg)$adj.r.squared
}

#sprat
idxF <- match(YEAR,sprat$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,spratLF$Year)
B <- spratLF[idxB+1,-1] #LF w/o year

year[26,]<-c(min(sprat$Year[idxF]),max(sprat$Year[idxF]))

size <- fishsize(B,seq(5,245,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[26,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[26,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[26,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[26,j] <- summary(Breg)$adj.r.squared
}

#whiting
idxF <- match(YEAR,whiting$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

idxB<- match(YEAR,spratLF$Year)
B <- whitingLF[idxB+1,-1] #LF w/o year

year[27,]<-c(min(whiting$Year[idxF]),max(whiting$Year[idxF]))

size <- fishsize(B,seq(10,690,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_1yr[27,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_1yr[27,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_1yr[27,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_1yr[27,j] <- summary(Breg)$adj.r.squared
}

colnames(sigF_1yr) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSST_1yr) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSSTF_1yr) <- c("mean","L95","shannon","evenness","skewness")

write.csv(sigF_1yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigF_1yr.csv")
write.csv(sigSST_1yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigSST_1yr.csv")
write.csv(sigSSTF_1yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigSSTF_1yr.csv")
write.csv(r2B_1yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_r2B_1yr.csv")
write.csv(year,file="~/Desktop/fishsize_varpart/output/year_analyzed.csv")

# 3-yr lag
colnames(sst)<-c("YEAR","SST")
sigSST_3yr <- matrix(NA,27,5)
sigF_3yr <- matrix(NA,27,5)
sigSSTF_3yr <- matrix(NA,27,5)
r2B_3yr <- matrix(NA,27,5)  
year <- matrix(NA,27,2)

# West US
# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[1,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(12,80,2))
L95 <- L95perc(B,size)
mean <- meansize(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[1,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[1,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[1,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[1,j] <- summary(Breg)$adj.r.squared
}

# Chilipepper RF
B <- chilipepperRFfnum
idxF <- match(chilipepperRFf$YEAR,exploitation$YEAR)
idxSST <- match(chilipepperRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[2,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(16,52,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[2,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[2,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[2,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[2,j] <- summary(Breg)$adj.r.squared
}

#Darkblotched RF
B <- darkblotchedRFfnum
idxF <- match(darkblotchedRFf$YEAR,exploitation$YEAR)
idxSST <- match(darkblotchedRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[3,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,c(6:33,seq(35,51,2)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[3,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[3,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[3,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[3,j] <- summary(Breg)$adj.r.squared
}

#Dover sole
B <- dovernum
idxF <- match(doverS$YEAR,exploitation$YEAR)
idxSST <- match(doverS$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dover.Sole[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[4,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(22,60,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[4,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[4,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[4,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[4,j] <- summary(Breg)$adj.r.squared
}

#English sole
B <- englishFnum
idxF <- match(englishF$YEAR,exploitation$YEAR)
idxSST <- match(englishF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$English.Sole[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[5,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(11,45,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[5,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[5,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[5,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[5,j] <- summary(Breg)$adj.r.squared
}

#Lingcod
B <- lingcodFnum
idxF <- match(lingcodF$YEAR,exploitation$YEAR)
idxSST <- match(lingcodF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Lingcod[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[6,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(28,110,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[6,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[6,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[6,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[6,j] <- summary(Breg)$adj.r.squared
}

#Longspine
B <- longspineNum
idxF <- match(longspine$YEAR,exploitation$YEAR)
idxSST <- match(longspine$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[7,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(5,35,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[7,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[7,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[7,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[7,j] <- summary(Breg)$adj.r.squared
}

#Petrale Sole
B <- petraleSoleFnum
idxF <- match(petraleSole$YEAR,exploitation$YEAR)
idxSST <- match(petraleSole$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Petrale.sole[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[8,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(12,62,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[8,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[8,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[8,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[8,j] <- summary(Breg)$adj.r.squared
}

#Sardine
B <- sardineNum[-4,] #exclude 1984
idxF <- match(sardineS1$YEAR,exploitation$YEAR)
idxSST <- match(sardineS1$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Sardine[idxF-3],sst$SST[idxSST-3]))[-4,]
colnames(A)<-c("F","SST")
year[9,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(9,26,0.5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[9,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[9,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[9,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[9,j] <- summary(Breg)$adj.r.squared
}

#Splitnose RF
B <- splitnoseRFfnum
idxF <- match(splitnoseRFf$YEAR,exploitation$YEAR)
idxSST <- match(splitnoseRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[10,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(2,56,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[10,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[10,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[10,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[10,j] <- summary(Breg)$adj.r.squared
}

#Yelloweye RF
B <- yelloweyeRFnum
idxF <- match(yelloweyeRF$YEAR,exploitation$YEAR)
idxSST <- match(yelloweyeRF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF-3],sst$SST[idxSST-3]))
colnames(A)<-c("F","SST")
year[11,] <- c(min(sst$YEAR[idxSST-3]),max(sst$YEAR[idxSST-3]))

size <- fishsize(B,seq(16,88,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[11,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[11,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[11,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[11,j] <- summary(Breg)$adj.r.squared
}

# Alaska
#AI pollock
B <- AIpollocknum
idxF <- match(AIpollockLF$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxSST <- match(exploitation$YEAR[idxF],AIsst$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..AI.[idxF-3],AIsst$SST[idxSST-3]))
idxB <- match(exploitation$YEAR[idxF],AIpollockLF$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[12,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(10,104,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[12,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[12,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[12,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[12,j] <- summary(Breg)$adj.r.squared
}

#GOA pollock
B <- GOApollockNum
idxF <- match(GOApollock$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..GOA.[idxF-3],GAK1bottomT$T[idxGAK1-3]))
idxB <- match(exploitation$YEAR[idxF],GOApollock$YEAR)
B <- B[idxB,]
year[13,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

colnames(A)<-c("F","SST")
size <- fishsize(B,seq(5,75,1))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[13,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[13,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[13,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[13,j] <- summary(Breg)$adj.r.squared
}

#BSAI flathead sole
#Check again
year<-seq(min(summerBottomT$YEAR),max(BSAIflatheadSole$YEAR)-3)
B <- BSAIflatheadSoleNum
idxsummerBottomT <- match(year,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(year,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Flathead.sole..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

year<-seq(min(summerBottomT$YEAR)+3,max(BSAIflatheadSole$YEAR))
idxB <- match(year,BSAIflatheadSoleF$YEAR)
B <- B[idxB,]  #shift the length distribution instead of temperature and F due to data availability
colnames(A)<-c("F","SST")
#year[14,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,c(seq(6,40,2),seq(43,58,3)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[14,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[14,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[14,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[14,j] <- summary(Breg)$adj.r.squared
}

#GOA flathead sole (DONE)
B <- GOAflatheadSoleNum
idxF <- match(GOAflatheadSoleF$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Flathead.sole..GOA.[idxF-3],GAK1bottomT$T[idxGAK1-3]))

idxB <- match(exploitation$YEAR[idxF],GOAflatheadSole$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
year[15,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(14,48,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[15,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[15,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[15,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[15,j] <- summary(Breg)$adj.r.squared
}

#BSAI pacific cod 
year<-seq(min(summerBottomT$YEAR),max(BSAIpacifcCod$YEAR)-3)
B <- BSAIpacifcCodNum
idxsummerBottomT <- match(year,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(year,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Pacific.cod..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

year<-seq(min(summerBottomT$YEAR)+3,max(BSAIpacifcCod$YEAR))
idxB <- match(year,BSAIpacifcCod$YEAR)
B <- B[idxB,]
colnames(A)<-c("F","SST")
#year[16,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,c(seq(9,42,3),seq(45,105,5)))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[16,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[16,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[16,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[16,j] <- summary(Breg)$adj.r.squared
}

#GOA rex sole
B <- GOArexSoleNum
idxF <- match(GOArexSoleF$YEAR,exploitation$YEAR)
idxGAK1 <- match(GOArexSoleF$YEAR,GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Rex.sole..GOA.[idxF-3],GAK1bottomT$T[idxGAK1-3]))
colnames(A)<-c("F","SST")
#year[17,] <- c(min(sst$YEAR[idxSST]),max(sst$YEAR[idxSST]))

size <- fishsize(B,seq(9,41,2))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$SST+A$F+A$SST*A$F)
  sigSST_3yr[17,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[17,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[17,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[17,j] <- summary(Breg)$adj.r.squared
}

# North Sea
#cod
YEAR<-1977:2011
idxF <- match(YEAR,cod$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

#year[18,]<-c(min(cod$Year[idxF]),max(cod$Year[idxF]))

size <- fishsize(B,seq(10,1400,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[18,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[18,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[18,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[18,j] <- summary(Breg)$adj.r.squared
}

#haddock
YEAR<-1977:2011
idxF <- match(YEAR,haddock$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

#year[19,]<-c(min(haddock$Year[idxF]),max(haddock$Year[idxF]))

size <- fishsize(B,seq(10,870,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[19,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[19,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[19,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[19,j] <- summary(Breg)$adj.r.squared
}

#herring
YEAR<-1977:2011
idxF <- match(YEAR,herring$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

#year[20,]<-c(min(herring$Year[idxF-1]),max(herring$Year[idxF-1]))

size <- fishsize(B,seq(5,380,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[20,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[20,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[20,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[20,j] <- summary(Breg)$adj.r.squared
}

#mackerel
YEAR<-1977:2011
idxF <- match(YEAR,mackerel$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(mackerel$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

#year[21,]<-c(min(mackerel$Year[idxF]),max(mackerel$Year[idxF]))

size <- fishsize(B,seq(10,580,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[21,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[21,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[21,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[21,j] <- summary(Breg)$adj.r.squared
}

#norwaypout
YEAR<-1977:2011
idxF <- match(YEAR,norwaypout$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(norwaypout$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

year[22,]<-c(min(norwaypout$Year[idxF]),max(norwaypout$Year[idxF]))

size <- fishsize(B,seq(10,320,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[22,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[22,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[22,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[22,j] <- summary(Breg)$adj.r.squared
}

#plaice
YEAR<-1977:2011
idxF <- match(YEAR,plaice$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

#year[23,]<-c(min(plaice$Year[idxF]),max(plaice$Year[idxF]))

size <- fishsize(B,seq(10,650,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[23,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[23,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[23,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[23,j] <- summary(Breg)$adj.r.squared
}

#saithe
YEAR<-1977:2011
idxF <- match(YEAR,saithe$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

#year[24,]<-c(min(saithe$Year[idxF]),max(saithe$Year[idxF]))

size <- fishsize(B,seq(10,1180,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[24,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[24,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[24,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[24,j] <- summary(Breg)$adj.r.squared
}

#sole
YEAR<-1977:2011
idxF <- match(YEAR,sole$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

year[25,]<-c(min(sole$Year[idxF]),max(sole$Year[idxF]))

size <- fishsize(B,seq(10,510,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[25,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[25,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[25,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[25,j] <- summary(Breg)$adj.r.squared
}

#sprat
YEAR<-1977:2011
idxF <- match(YEAR,sprat$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")

YEAR<-1980:2014
idxB<- match(YEAR,spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

#year[26,]<-c(min(sprat$Year[idxF]),max(sprat$Year[idxF]))

size <- fishsize(B,seq(5,245,5))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[26,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[26,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[26,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[26,j] <- summary(Breg)$adj.r.squared
}

#whiting
YEAR<-1977:2011
idxF <- match(YEAR,whiting$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("F","bottomT")
YEAR<-1980:2014
idxB<- match(YEAR,spratLF$Year)
B <- whitingLF[idxB,-1] #LF w/o year

#year[27,]<-c(min(whiting$Year[idxF]),max(whiting$Year[idxF]))

size <- fishsize(B,seq(10,690,10))
mean <- meansize(B,size)
L95 <- L95perc(B,size)
shan <- shannon(B)
even <- evenness(B)
skew <- skewness(B)

newB <- data.frame(cbind(mean,L95,shan,even,skew))

for (j in 1:ncol(newB)){
  Breg <- lm(newB[,j]~1+A$bottomT+A$F+A$bottomT*A$F)
  sigSST_3yr[27,j] <- anova(Breg)$`Pr(>F)`[1]
  sigF_3yr[27,j] <- anova(Breg)$`Pr(>F)`[2]
  sigSSTF_3yr[27,j] <- anova(Breg)$`Pr(>F)`[3]
  r2B_3yr[27,j] <- summary(Breg)$adj.r.squared
}

colnames(sigF_3yr) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSST_3yr) <- c("mean","L95","shannon","evenness","skewness")
colnames(sigSSTF_3yr) <- c("mean","L95","shannon","evenness","skewness")

write.csv(sigF_3yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigF_3yr.csv")
write.csv(sigSST_3yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigSST_3yr.csv")
write.csv(sigSSTF_3yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_sigSSTF_3yr.csv")
write.csv(r2B_3yr,file = "~/Desktop/fishsize_varpart/output/uniSBI_r2B_3yr.csv")
write.csv(year,file="~/Desktop/fishsize_varpart/output/year_analyzed.csv")