# Variance partitioning
# 20160929 Crystal 
# 20180211 Use survey length, bottom temperature and transformed fishing mortality
# 20180225 Incl. timelag
library(vegan)
library(spdep)

meanA <- matrix(NA,28,2) # West US: 12 species, Alaska: 7 stocks, North Sea: 9 species
cvA <- matrix(NA,28,2)

lifehist <- read.csv("~/Desktop/fishsize_varpart/data/lifehist_habitat.csv")

# West US
sst <- read.table("~/Desktop/fishsize_varpart/data/WestUS/WestUSannualSSTfrom1948.txt",header = FALSE)
exploitation <- read.csv("~/Desktop/fishsize_varpart/data/FishingMortality.csv",header=TRUE)
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

#AI pollock
ai_pollockAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/ai_pollockAnnual.csv",header=TRUE)
AIbottomT <- ai_pollockAnnual[,c(1,5)]
  
B <- AIpollocknum
idxF <- match(AIbottomT$YEAR,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxTEMP <- match(AIbottomT$YEAR,exploitation$YEAR[idxF])
idxTEMP <- idxTEMP[is.finite(idxTEMP)]

A <- data.frame(cbind(exploitation$Walleye.pollock..AI.[idxF],AIbottomT$BOT_TEMP[idxTEMP]))
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
goa_pollockAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/goa_pollockAnnual.csv",header=TRUE)
B <- GOApollockNum

goaBottomT <- goa_pollockAnnual[,c(1,5)]
idxBottomT <- match(GOApollock$YEAR,goaBottomT$YEAR)
idxBottomT <- idxBottomT[is.finite(idxBottomT)]

idxF <- match(goaBottomT$YEAR[idxBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxB <- match(goaBottomT$YEAR[idxBottomT],GOApollock$YEAR)
idxB <- idxB[is.finite(idxB)]
B <- B[idxB,]
A <- data.frame(cbind(exploitation$Walleye.pollock..GOA.[idxF],goaBottomT$BOT_TEMP[idxBottomT]))

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

#EBS flathead sole
ebs_flatheadSoleAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/ebs_flatheadSoleAnnual.csv",header = TRUE)
ebsBottomT <- ebs_flatheadSoleAnnual[,c(1,5)]

B <- BSAIflatheadSoleNum
idxBottomT <- match(BSAIflatheadSoleF$YEAR,ebsBottomT$YEAR)
idxBottomT <- idxBottomT[is.finite(idxBottomT)]
idxF <- match(ebsBottomT$YEAR[idxBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Flathead.sole.BSAI[idxF],ebsBottomT$BOT_TEMP[idxBottomT-1]))

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
goa_flatheadSoleAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/goa_flatheadSoleAnnual.csv",header=TRUE)
goaBottomT <- goa_flatheadSoleAnnual[,c(1,5)]

B <- GOAflatheadSoleNum
idxF <- match(GOAflatheadSoleF$Year,exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
idxBottomT <- match(exploitation$YEAR[idxF],goaBottomT$YEAR)
A <- data.frame(cbind(exploitation$Flathead.sole.GOA[idxF],goaBottomT$BOT_TEMP[idxBottomT]))

idxB <- match(exploitation$YEAR[idxF],GOAflatheadSole$Year)
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

#EBS pacific cod
ebs_pacificCodAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/ebs_pacificCodAnnual.csv",header = TRUE)
ebsBottomT <- ebs_pacificCodAnnual[,c(1,5)]

B <- EBSpacifcCodNum

idxBottomT <- match(EBSpacifcCod$YEAR,ebsBottomT$YEAR)
idxBottomT <- idxBottomT[is.finite(idxBottomT)]
idxF <- match(ebsBottomT$YEAR[idxBottomT],exploitation$YEAR)
idxF <- idxF[is.finite(idxF)]
A <- data.frame(cbind(exploitation$Pacific.cod.BSAI[idxF],ebsBottomT$BOT_TEMP[idxBottomT]))

idxB <- match(exploitation$YEAR[idxF],EBSpacifcCod$YEAR)
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

#GOA pacific cod
goa_pacificCodAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/goa_pacificCodAnnual.csv",header=TRUE)
goaBottomT <- goa_pacificCodAnnual[,c(1,5)]

B <- GOApacificCodNum

idxF <- match(GOApacificCod$YEAR,exploitation$YEAR)
idxBottomT <- match(GOApacificCod$YEAR,goaBottomT$YEAR)
A <- data.frame(cbind(exploitation$Pacific.cod.GOA[idxF],goaBottomT$BOT_TEMP[idxBottomT]))

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

#GOA rex sole
goa_rexSoleAnnual <- read.csv("~/Desktop/fishsize_varpart/data/Alaska/envi/goa_rexSoleAnnual.csv",header=TRUE)
goaBottomT <- goa_rexSoleAnnual[,c(1,5)]

B <- GOArexSoleNum

idxF <- match(GOArexSoleF$YEAR,exploitation$YEAR)
idxBottomT <- match(GOArexSoleF$YEAR,goaBottomT$YEAR)
A <- data.frame(cbind(exploitation$Rex.sole.GOA[idxF],goaBottomT$BOT_TEMP[idxBottomT]))

M <- lifehist$M[18]

meanA[18,] <- c(mean(A$X1/M),mean(A$X2))
cvA[18,] <- c(sd(A$X1)/mean(A$X1),sd(A$X2)/mean(A$X2))

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
NSavg_tempQ1 <- read.csv("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/NSavg_tempQ1.csv",header = TRUE)
year<-1977:2014

#cod
idxF <- match(year,cod$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

M <- cod$M[idxF]

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

#haddock
idxF <- match(year,haddock$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

M <- haddock$M[idxF]

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

#herring
idxF <- match(year,herring$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

M <- herring$M[idxF]

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

#mackerel 1980-
idxF <- match(year,mackerel$Year)
idxbottomT <- match(mackerel$Year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(mackerel$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

idxB<- match(mackerel$Year[idxF[!is.na(idxF)]],mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

M <- lifehist$M[22]

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

#norwaypout 1984-
idxF <- match(year,norwaypout$Year)
idxbottomT <- match(norwaypout$Year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(norwaypout$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

idxB<- match(norwaypout$Year[idxF[!is.na(idxF)]],norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

M <- norwaypout$M[idxF[!is.na(idxF)]]

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

#plaice
idxF <- match(year,plaice$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

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

#saithe
idxF <- match(year,saithe$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

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

#sole
idxF <- match(year,sole$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(year,soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

M <- lifehist$M[26]

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

#sprat
idxF <- match(year,sprat$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

M <- sprat$M[idxF]

idxB<- match(year,spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

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

#whiting 1990-
idxF <- match(year,whiting$Year)
idxbottomT <- match(whiting$Year,NSavg_tempQ1$YEAR)
A <- data.frame(cbind(whiting$F[idxF[!is.na(idxF)]],NSavg_tempQ1$temperature[idxbottomT[!is.na(idxbottomT)]]))
colnames(A)=c("exploitation","temperature")

M <- whiting$M[idxF[!is.na(idxF)]]

idxB<- match(whiting$Year[idxF[!is.na(idxF)]],whitingLF$Year)
B <- whitingLF[idxB,-1] #LF w/o year

meanA[28,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA[28,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

adjR2Bout[adjR2Bout<0]=0

write.csv(adjR2Bout,file = "output/adjR2Bout.csv")
write.csv(sigout,file = "output/sigout.csv")
write.csv(meanA,file="output/meanSST_F.csv")
write.csv(cvA,file="output/cvSST_F.csv")

rm(adjR2Bout,sigout)