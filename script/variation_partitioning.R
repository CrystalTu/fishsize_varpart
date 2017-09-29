# Variance partitioning
# 20160929 Crystal 

library(vegan)

meanA <- matrix(NA,27,2) # West US: 12 species, Alaska: 6 species, North Sea: 9 species
cvA <- matrix(NA,27,2)

lifehist <- read.csv("lifehist_habitat.csv")

# West US
source("~/Desktop/fishsize_varpart/data/length/preset.R") #length 
sst <- read.table("~/Desktop/fishsize_varpart/data/WestUS/WestUSannualSSTfrom1948.txt",header = FALSE)
exploitation <- read.csv("~/Desktop/fishsize_varpart/data/WestUS/exploitation.csv",header=TRUE)
colnames(sst)<-c("YEAR","SST")

# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF],sst$SST[idxSST]))

M <- lifehist$M[1]

meanA[1,] <- c(mean(A$X1/M),mean(A$X2))
cvA[1,] <- c(sd(A$X1/M)/mean(A$X1/M),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

# output
adjR2Bout <- adjR2B
sigout <- sig

# Chilipepper RF
B <- chilipepperRFfnum
idxF <- match(chilipepperRFf$YEAR,exploitation$YEAR)
idxSST <- match(chilipepperRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF],sst$SST[idxSST]))

M <- lifehist$M[2]

meanA[2,] <- c(mean(A$X1/M),mean(A$X2))
cvA[2,] <- c(sd(A$X1/lifehist$M[2])/mean(A$X1/lifehist$M[2]),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Darkblotched RF
B <- darkblotchedRFfnum
idxF <- match(darkblotchedRFf$YEAR,exploitation$YEAR)
idxSST <- match(darkblotchedRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF],sst$SST[idxSST]))

M <- lifehist$M[3]

meanA[3,] <- c(mean(A$X1/M),mean(A$X2))
cvA[3,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Dover sole
B <- dovernum
idxF <- match(doverS$YEAR,exploitation$YEAR)
idxSST <- match(doverS$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dover.Sole[idxF],sst$SST[idxSST]))

M <- lifehist$M[4]

meanA[4,] <- c(mean(A$X1/M),mean(A$X2))
cvA[4,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#English sole
B <- englishFnum
idxF <- match(englishF$YEAR,exploitation$YEAR)
idxSST <- match(englishF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$English.Sole[idxF],sst$SST[idxSST]))

M <- lifehist$M[5]

meanA[5,] <- c(mean(A$X1/M),mean(A$X2))
cvA[5,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Lingcod
B <- lingcodFnum
idxF <- match(lingcodF$YEAR,exploitation$YEAR)
idxSST <- match(lingcodF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Lingcod[idxF],sst$SST[idxSST]))

M <- lifehist$M[6]

meanA[6,] <- c(mean(A$X1/M),mean(A$X2))
cvA[6,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Longspine
B <- longspineNum
idxF <- match(longspine$YEAR,exploitation$YEAR)
idxSST <- match(longspine$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF],sst$SST[idxSST]))

M <- lifehist$M[7]

meanA[7,] <- c(mean(A$X1/M),mean(A$X2))
cvA[7,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Petrale Sole
B <- petraleSoleFnum
idxF <- match(petraleSole$YEAR,exploitation$YEAR)
idxSST <- match(petraleSole$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Petrale.sole[idxF],sst$SST[idxSST]))

M <- lifehist$M[8]

meanA[8,] <- c(mean(A$X1/M),mean(A$X2))
cvA[8,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Sardine
B <- sardineNum[-4,] #exclude 1984
idxF <- match(sardineS1$YEAR,exploitation$YEAR)
idxSST <- match(sardineS1$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Sardine[idxF],sst$SST[idxSST]))[-4,]

M <- lifehist$M[9]

meanA[9,] <- c(mean(A$X1/M),mean(A$X2))
cvA[9,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A) 
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Splitnose RF
B <- splitnoseRFfnum
idxF <- match(splitnoseRFf$YEAR,exploitation$YEAR)
idxSST <- match(splitnoseRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF],sst$SST[idxSST]))

M <- lifehist$M[10]

meanA[10,] <- c(mean(A$X1/M),mean(A$X2))
cvA[10,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Yelloweye RF
B <- yelloweyeRFnum
idxF <- match(yelloweyeRF$YEAR,exploitation$YEAR)
idxSST <- match(yelloweyeRF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF],sst$SST[idxSST]))

M <- lifehist$M[11]

meanA[11,] <- c(mean(A$X1/M),mean(A$X2))
cvA[11,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

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

M <- lifehist$M[12]

meanA[12,] <- c(mean(A$X1/M),mean(A$X2))
cvA[12,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA pollock
B <- GOApollockNum
idxF <- match(GOApollock$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..GOA.[idxF],GAK1bottomT$T[idxGAK1]))
idxB <- match(exploitation$YEAR[idxF],GOApollock$YEAR)
B <- B[idxB,]

M <- lifehist$M[13]

meanA[13,] <- c(mean(A$X1/M),mean(A$X2))
cvA[13,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#BSAI flathead sole (done)
B <- BSAIflatheadSoleNum
idxsummerBottomT <- match(BSAIflatheadSoleF$YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Flathead.sole..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT-1]))

idxB <- match(exploitation$YEAR[idxF],BSAIflatheadSoleF$YEAR)
B <- B[idxB,]

M <- lifehist$M[14]

meanA[14,] <- c(mean(A$X1/M),mean(A$X2))
cvA[14,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA flathead sole (DONE)
B <- GOAflatheadSoleNum
idxF <- match(GOAflatheadSole$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Flathead.sole..GOA.[idxF],GAK1bottomT$T[idxGAK1]))

idxB <- match(exploitation$YEAR[idxF],GOAflatheadSole$YEAR)
B <- B[idxB,]

M <- lifehist$M[15]

meanA[15,] <- c(mean(A$X1/M),mean(A$X2))
cvA[15,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#BSAI pacific cod
B <- BSAIpacifcCodNum
idxsummerBottomT <- match(BSAIpacifcCod$YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Pacific.cod..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

idxB <- match(exploitation$YEAR[idxF],BSAIpacifcCod$YEAR)
B <- B[idxB,]

M <- lifehist$M[16]

meanA[16,] <- c(mean(A$X1/M),mean(A$X2))
cvA[16,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA rex sole
B <- GOArexSoleNum
idxF <- match(GOArexSoleF$YEAR,exploitation$YEAR)
idxGAK1 <- match(GOArexSoleF$YEAR,GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Rex.sole..GOA.[idxF],GAK1bottomT$T[idxGAK1]))

M <- lifehist$M[17]

meanA[17,] <- c(mean(A$X1/M),mean(A$X2))
cvA[17,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

# North Sea
setwd("~/Desktop/fishsize_varpart/")
source("~/Desktop/fishsize_varpart/data/NorthSea/IBTSsurvey/preset.R") # get length frequency data
source("~/Desktop/fishsize_varpart/data/NorthSea/stockassessment/stockassessment.R") # get stockassessment result table
source("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/temperature.R") # get bottom T
setwd("~/Desktop/fishsize_varpart/")

year<-1977:2014

#cod
idxF <- match(year,cod$Year)
idxbottomT <- match(year,NSavg_tempQ1$Year)
A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

M <- cod$M[idxF]

meanA[18,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[18,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#haddock
idxF <- match(year,haddock$Year)
idxbottomT <- match(year,NSavg_tempQ1$Year)
A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

M <- haddock$M[idxF]

meanA[19,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[19,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#herring
idxF <- match(year,herring$Year)
idxbottomT <- match(year,NSavg_tempQ1$Year)
A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

M <- herring$M[idxF]

meanA[20,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[20,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#mackerel 1980-
idxF <- match(year,mackerel$Year)
idxbottomT <- match(mackerel$Year,NSavg_tempQ1$Year)
A <- data.frame(cbind(mackerel$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

idxB<- match(mackerel$Year[idxF[!is.na(idxF)]],mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

M <- lifehist$M[21]

meanA[21,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[21,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#norwaypout 1984-
idxF <- match(year,norwaypout$Year)
idxbottomT <- match(norwaypout$Year,NSavg_tempQ1$Year)
A <- data.frame(cbind(norwaypout$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

idxB<- match(norwaypout$Year[idxF[!is.na(idxF)]],norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

M <- norwaypout$M[idxF[!is.na(idxF)]]

meanA[22,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[22,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#plaice
idxF <- match(year,plaice$Year)
idxbottomT <- match(year,NSavg_tempQ1$Year)
A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

M <- lifehist$M[23]

meanA[23,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[23,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#saithe
idxF <- match(year,saithe$Year)
idxbottomT <- match(year,NSavg_tempQ1$Year)
A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

M <- lifehist$M[24]

meanA[24,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[24,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#sole
idxF <- match(year,sole$Year)
idxbottomT <- match(year,NSavg_tempQ1$Year)
A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

M <- lifehist$M[25]

meanA[25,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[25,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#sprat
idxF <- match(year,sprat$Year)
idxbottomT <- match(year,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

M <- sprat$M[idxF]

idxB<- match(year,spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

meanA[26,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[26,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#whiting 1990-
idxF <- match(year,whiting$Year)
idxbottomT <- match(whiting$Year,NSavg_tempQ1$Year)
A <- data.frame(cbind(whiting$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

M <- whiting$M[idxF[!is.na(idxF)]]

idxB<- match(whiting$Year[idxF[!is.na(idxF)]],whitingLF$Year)
B <- whitingLF[idxB,-1] #LF w/o year

meanA[27,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[27,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

adjR2Bout[adjR2Bout<0]=0

write.csv(adjR2Bout,file = "~/Desktop/fishsize_varpart/output/adjR2Bout.csv")
write.csv(sigout,file = "~/Desktop/fishsize_varpart/output/sigout.csv")
write.csv(meanA,file="~/Desktop/fishsize_varpart/output/meanSST_F.csv")
write.csv(cvA,file="~/Desktop/fishsize_varpart/output/cvSST_F.csv")

rm(adjR2Bout,sigout)

# 1-year lag

meanA1yr <- matrix(NA,27,2) # West US: 12 species, Alaska: 6 species, North Sea: 9 species
cvA1yr <- matrix(NA,27,2)

lag=1

library(vegan)
# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[1,] <- colMeans(A)
cvA1yr[1,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

# output
adjR2Bout <- adjR2B
sigout <- sig

# Chilipepper RF
B <- chilipepperRFfnum
idxF <- match(chilipepperRFf$YEAR,exploitation$YEAR)
idxSST <- match(chilipepperRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[2,] <- colMeans(A)
cvA1yr[2,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Darkblotched RF
B <- darkblotchedRFfnum
idxF <- match(darkblotchedRFf$YEAR,exploitation$YEAR)
idxSST <- match(darkblotchedRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[3,] <- colMeans(A)
cvA1yr[3,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Dover sole
B <- dovernum
idxF <- match(doverS$YEAR,exploitation$YEAR)
idxSST <- match(doverS$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dover.Sole[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[4,] <- colMeans(A)
cvA1yr[4,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#English sole
B <- englishFnum
idxF <- match(englishF$YEAR,exploitation$YEAR)
idxSST <- match(englishF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$English.Sole[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[5,] <- colMeans(A)
cvA1yr[5,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Lingcod
B <- lingcodFnum
idxF <- match(lingcodF$YEAR,exploitation$YEAR)
idxSST <- match(lingcodF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Lingcod[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[6,] <- colMeans(A)
cvA1yr[6,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Longspine
B <- longspineNum
idxF <- match(longspine$YEAR,exploitation$YEAR)
idxSST <- match(longspine$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[7,] <- colMeans(A)
cvA1yr[7,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Petrale Sole
B <- petraleSoleFnum
idxF <- match(petraleSole$YEAR,exploitation$YEAR)
idxSST <- match(petraleSole$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Petrale.sole[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[8,] <- colMeans(A)
cvA1yr[8,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Sardine
B <- sardineNum[-4,] #exclude 1984
idxF <- match(sardineS1$YEAR[-4],exploitation$YEAR)
idxSST <- match(sardineS1$YEAR[-4],sst$YEAR)
A <- data.frame(cbind(exploitation$Sardine[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[9,] <- colMeans(A)
cvA1yr[9,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

A <- A[-1,]
B <- B[-1,] #match the dimension

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A) 
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Splitnose RF
B <- splitnoseRFfnum
idxF <- match(splitnoseRFf$YEAR,exploitation$YEAR)
idxSST <- match(splitnoseRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[10,] <- colMeans(A)
cvA1yr[10,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Yelloweye RF
B <- yelloweyeRFnum
idxF <- match(yelloweyeRF$YEAR,exploitation$YEAR)
idxSST <- match(yelloweyeRF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA1yr[11,] <- colMeans(A)
cvA1yr[11,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Alaska
#AI pollock
B <- AIpollocknum
idxF <- match(AIpollockLF$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxSST <- match(exploitation$YEAR[idxF],AIsst$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..AI.[idxF-lag],AIsst$SST[idxSST-lag]))
idxB <- match(exploitation$YEAR[idxF],AIpollockLF$YEAR)
B <- B[idxB,]

A <- A[-1,]
B <- B[-1,] #match the dimension

meanA1yr[12,] <- colMeans(A)
cvA1yr[12,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA pollock
B <- GOApollockNum
idxF <- match(GOApollock$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..GOA.[idxF-lag],GAK1bottomT$T[idxGAK1-lag]))
idxB <- match(exploitation$YEAR[idxF],GOApollock$YEAR)
B <- B[idxB,]
meanA1yr[13,] <- colMeans(A)
cvA1yr[13,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#BSAI flathead sole (done)
YEAR<-seq(min(summerBottomT$YEAR),max(BSAIpacifcCod$YEAR)-1)
B <- BSAIflatheadSoleNum
idxsummerBottomT <- match(YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Flathead.sole..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

YEAR<-seq(min(summerBottomT$YEAR)+1,max(BSAIpacifcCod$YEAR))
idxB <- match(YEAR,BSAIpacifcCod$YEAR)
B <- B[idxB,]

meanA1yr[14,] <- colMeans(A)
cvA1yr[14,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA flathead sole (DONE)
B <- GOAflatheadSoleNum
idxF <- match(GOAflatheadSole$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Flathead.sole..GOA.[idxF-lag],GAK1bottomT$T[idxGAK1-lag]))

idxB <- match(exploitation$YEAR[idxF],GOAflatheadSole$YEAR)
B <- B[idxB,]

A <- A[-c(1,27),]
B <- B[-c(1,27),] #match the dimension

meanA1yr[15,] <- colMeans(A)
cvA1yr[15,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

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
meanA1yr[16,] <- colMeans(A)
cvA1yr[16,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA rex sole
B <- GOArexSoleNum
idxF <- match(GOArexSoleF$YEAR,exploitation$YEAR)
idxGAK1 <- match(GOArexSoleF$YEAR,GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Rex.sole..GOA.[idxF-lag],GAK1bottomT$T[idxGAK1-lag]))
A <- A[-1,] #exploitation for 1981 is not available
B <- B[-1,] #match accordingly (1983-)

meanA1yr[17,] <- colMeans(A)
cvA1yr[17,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#North Sea 20170205
year<-1977:2014

#cod
YEAR<-1977:2013
idxF <- match(YEAR,cod$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1978:2014
idxB<- match(YEAR,codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

meanA1yr[18,] <- colMeans(A)
cvA1yr[18,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#haddock
YEAR<-1977:2013
idxF <- match(YEAR,haddock$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1978:2014
idxB<- match(YEAR,haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

meanA1yr[19,] <- colMeans(A)
cvA1yr[19,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#herring
YEAR<-1977:2013
idxF <- match(YEAR,herring$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1978:2014
idxB<- match(YEAR,herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

meanA1yr[20,] <- colMeans(A)
cvA1yr[20,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#mackerel 1980-
YEAR<-1980:2013
idxF <- match(YEAR,mackerel$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(mackerel$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

YEAR<-1981:2014
idxB<- match(YEAR,mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

meanA1yr[21,] <- colMeans(A)
cvA1yr[21,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#norwaypout 1984-
YEAR<-1984:2013
idxF <- match(YEAR,norwaypout$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(norwaypout$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

YEAR<-1985:2014
idxB<- match(YEAR,norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

meanA1yr[22,] <- colMeans(A)
cvA1yr[22,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#plaice
YEAR<-1977:2013
idxF <- match(YEAR,plaice$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1978:2014
idxB<- match(YEAR,plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

meanA1yr[23,] <- colMeans(A)
cvA1yr[23,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#saithe
YEAR<-1977:2013
idxF <- match(YEAR,saithe$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1978:2014
idxB<- match(YEAR,saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

meanA1yr[24,] <- colMeans(A)
cvA1yr[24,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#sole
YEAR<-1977:2013
idxF <- match(YEAR,sole$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1978:2014
idxB<- match(YEAR,soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

meanA1yr[25,] <- colMeans(A)
cvA1yr[25,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#sprat
YEAR<-1977:2013
idxF <- match(YEAR,sprat$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1978:2014
idxB<- match(YEAR,spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

meanA1yr[26,] <- colMeans(A)
cvA1yr[26,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#whiting 1990-
YEAR<-1990:2013
idxF <- match(YEAR,whiting$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(whiting$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

YEAR<-1991:2014
idxB<- match(YEAR,whitingLF$Year)
B <- whitingLF[idxB,-1] #LF w/o year

meanA1yr[27,] <- colMeans(A)
cvA1yr[27,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

adjR2Bout[adjR2Bout<0]=0

write.csv(adjR2Bout,file = "~/Desktop/fishsize_varpart/output/adjR2Bout1yr.csv")
write.csv(sigout,file = "~/Desktop/fishsize_varpart/output/sigout1yr.csv")
write.csv(meanA1yr,file="~/Desktop/fishsize_varpart/output/meanSST_F1yr.csv")
write.csv(cvA1yr,file="~/Desktop/fishsize_varpart/output/cvSST_F1yr.csv")

rm(adjR2Bout,sigout)

# 3-yr lag
meanA3yr <- matrix(NA,27,2) # West US: 12 species, Alaska: 6 species, North Sea: 9 species
cvA3yr <- matrix(NA,27,2)

lag=3

library(vegan)
# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[1,] <- colMeans(A)
cvA3yr[1,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

# output
adjR2Bout <- adjR2B
sigout <- sig

# Chilipepper RF
B <- chilipepperRFfnum
idxF <- match(chilipepperRFf$YEAR,exploitation$YEAR)
idxSST <- match(chilipepperRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[2,] <- colMeans(A)
cvA3yr[2,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Darkblotched RF
B <- darkblotchedRFfnum
idxF <- match(darkblotchedRFf$YEAR,exploitation$YEAR)
idxSST <- match(darkblotchedRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[3,] <- colMeans(A)
cvA3yr[3,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Dover sole
B <- dovernum
idxF <- match(doverS$YEAR,exploitation$YEAR)
idxSST <- match(doverS$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Dover.Sole[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[4,] <- colMeans(A)
cvA3yr[4,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#English sole
B <- englishFnum
idxF <- match(englishF$YEAR,exploitation$YEAR)
idxSST <- match(englishF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$English.Sole[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[5,] <- colMeans(A)
cvA3yr[5,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Lingcod
B <- lingcodFnum
idxF <- match(lingcodF$YEAR,exploitation$YEAR)
idxSST <- match(lingcodF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Lingcod[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[6,] <- colMeans(A)
cvA3yr[6,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Longspine
B <- longspineNum
idxF <- match(longspine$YEAR,exploitation$YEAR)
idxSST <- match(longspine$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[7,] <- colMeans(A)
cvA3yr[7,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Petrale Sole
B <- petraleSoleFnum
idxF <- match(petraleSole$YEAR,exploitation$YEAR)
idxSST <- match(petraleSole$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Petrale.sole[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[8,] <- colMeans(A)
cvA3yr[8,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Sardine
B <- sardineNum[-4,] #exclude 1984
idxF <- match(sardineS1$YEAR[-4],exploitation$YEAR)
idxSST <- match(sardineS1$YEAR[-4],sst$YEAR)
A <- data.frame(cbind(exploitation$Sardine[idxF-lag],sst$SST[idxSST-lag]))

A <- A[-c(1,2,3),]
B <- B[-c(1,2,3),] #match the dimension

meanA3yr[9,] <- colMeans(A)
cvA3yr[9,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))
colnames(A)=c("exploitation","temperature")

partB<-varpart(B, X= ~exploitation,~temperature,data=A) 
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Splitnose RF
B <- splitnoseRFfnum
idxF <- match(splitnoseRFf$YEAR,exploitation$YEAR)
idxSST <- match(splitnoseRFf$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[10,] <- colMeans(A)
cvA3yr[10,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Yelloweye RF
B <- yelloweyeRFnum
idxF <- match(yelloweyeRF$YEAR,exploitation$YEAR)
idxSST <- match(yelloweyeRF$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF-lag],sst$SST[idxSST-lag]))
meanA3yr[11,] <- colMeans(A)
cvA3yr[11,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#Alaska
source("~/Desktop/fishsize_varpart/data/Alaska/preset.R")
setwd("~/Desktop/fishsize_varpart/")

#AI pollock
B <- AIpollocknum
idxF <- match(AIpollockLF$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxSST <- match(exploitation$YEAR[idxF],AIsst$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..AI.[idxF-lag],AIsst$SST[idxSST-lag]))
idxB <- match(exploitation$YEAR[idxF],AIpollockLF$YEAR)
B <- B[idxB,]

A <- A[-c(1,2,3),]
B <- B[-c(1,2,3),] #match the dimension

meanA3yr[12,] <- colMeans(A)
cvA3yr[12,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA pollock
B <- GOApollockNum
idxF <- match(GOApollock$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Walleye.pollock..GOA.[idxF-lag],GAK1bottomT$T[idxGAK1-lag]))
idxB <- match(exploitation$YEAR[idxF],GOApollock$YEAR)
B <- B[idxB,]
meanA3yr[13,] <- colMeans(A)
cvA3yr[13,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#BSAI flathead sole (done)
YEAR<-seq(min(summerBottomT$YEAR),max(BSAIpacifcCod$YEAR)-lag)
B <- BSAIflatheadSoleNum
idxsummerBottomT <- match(YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Flathead.sole..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

YEAR<-seq(min(summerBottomT$YEAR)+lag,max(BSAIpacifcCod$YEAR))
idxB <- match(YEAR,BSAIpacifcCod$YEAR)
B <- B[idxB,]
meanA3yr[14,] <- colMeans(A)
cvA3yr[14,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA flathead sole (DONE)
B <- GOAflatheadSoleNum
idxF <- match(GOAflatheadSoleF$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxGAK1 <- match(exploitation$YEAR[idxF],GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Flathead.sole..GOA.[idxF-lag],GAK1bottomT$T[idxGAK1-lag]))

idxB <- match(exploitation$YEAR[idxF],GOAflatheadSole$YEAR)
B <- B[idxB,]

A <- A[-c(1,2,3),] #remove the NA year
B <- B[-c(1,2,3),]

meanA3yr[15,] <- colMeans(A)
cvA3yr[15,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#BSAI pacific cod
YEAR<-seq(min(summerBottomT$YEAR),max(BSAIpacifcCod$YEAR)-lag)

B <- BSAIpacifcCodNum
idxsummerBottomT <- match(YEAR,summerBottomT$YEAR)
idxsummerBottomT <- idxsummerBottomT[is.finite(idxsummerBottomT)]
idxF <- match(summerBottomT$YEAR[idxsummerBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Pacific.cod..BSAI.[idxF],summerBottomT$summerBottomT[idxsummerBottomT]))

YEAR<-seq(min(summerBottomT$YEAR)+lag,max(BSAIpacifcCod$YEAR))
idxB <- match(YEAR,BSAIpacifcCod$YEAR)
B <- B[idxB,]
meanA3yr[16,] <- colMeans(A)
cvA3yr[16,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#GOA rex sole
B <- GOArexSoleNum
idxF <- match(GOArexSoleF$YEAR,exploitation$YEAR)
idxGAK1 <- match(GOArexSoleF$YEAR,GAK1bottomT$YEAR)
A <- data.frame(cbind(exploitation$Rex.sole..GOA.[idxF-lag],GAK1bottomT$T[idxGAK1-lag]))
A <- A[-c(1,2,3),] 
B <- B[-c(1,2,3),] #match accordingly

meanA3yr[17,] <- colMeans(A)
cvA3yr[17,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

colnames(A)=c("exploitation","temperature")
partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#North Sea
year<-1977:2014

#cod
YEAR<-1977:2011
idxF <- match(YEAR,cod$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1980:2014
idxB<- match(YEAR,codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

meanA3yr[18,] <- colMeans(A)
cvA3yr[18,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#haddock
YEAR<-1977:2011
idxF <- match(YEAR,haddock$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1980:2014
idxB<- match(YEAR,haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

meanA3yr[19,] <- colMeans(A)
cvA3yr[19,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#herring
YEAR<-1977:2011
idxF <- match(YEAR,herring$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1980:2014
idxB<- match(YEAR,herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

meanA3yr[20,] <- colMeans(A)
cvA3yr[20,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#mackerel 1980-
YEAR<-1980:2013
idxF <- match(YEAR,mackerel$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(mackerel$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

YEAR<-1981:2014
idxB<- match(YEAR,mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

meanA3yr[21,] <- colMeans(A)
cvA3yr[21,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#norwaypout 1984-
YEAR<-1984:2013
idxF <- match(YEAR,norwaypout$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(norwaypout$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

YEAR<-1985:2014
idxB<- match(YEAR,norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

meanA3yr[22,] <- colMeans(A)
cvA3yr[22,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#plaice
YEAR<-1977:2011
idxF <- match(YEAR,plaice$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1980:2014
idxB<- match(YEAR,plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

meanA3yr[23,] <- colMeans(A)
cvA3yr[23,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#saithe
YEAR<-1977:2011
idxF <- match(YEAR,saithe$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1980:2014
idxB<- match(YEAR,saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

meanA3yr[24,] <- colMeans(A)
cvA3yr[24,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#sole
YEAR<-1977:2011
idxF <- match(YEAR,sole$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1980:2014
idxB<- match(YEAR,soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

meanA3yr[25,] <- colMeans(A)
cvA3yr[25,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#sprat
YEAR<-1977:2011
idxF <- match(YEAR,sprat$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

YEAR<-1980:2014
idxB<- match(YEAR,spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

meanA3yr[26,] <- colMeans(A)
cvA3yr[26,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

#whiting 1990-
YEAR<-1990:2011
idxF <- match(YEAR,whiting$Year)
idxbottomT <- match(YEAR,NSavg_tempQ1$Year)
A <- data.frame(cbind(whiting$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

YEAR<-1993:2014
idxB<- match(YEAR,whitingLF$Year)
B <- whitingLF[idxB,-1] #LF w/o year

meanA3yr[27,] <- colMeans(A)
cvA3yr[27,] <- c(sd(A$exploitation)/mean(A$exploitation),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

adjR2Bout[adjR2Bout<0]=0

write.csv(adjR2Bout,file = "~/Desktop/fishsize_varpart/output/adjR2Bout3yr.csv")
write.csv(sigout,file = "~/Desktop/fishsize_varpart/output/sigout3yr.csv")
write.csv(meanA1yr,file="~/Desktop/fishsize_varpart/output/meanSST_F3yr.csv")
write.csv(cvA1yr,file="~/Desktop/fishsize_varpart/output/cvSST_F3yr.csv")
