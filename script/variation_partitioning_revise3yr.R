# 3yr (time-interval) lag for temperature (except AI and GOA)
# 20180225

meanA3yr <- matrix(NA,28,2) # West US: 12 species, Alaska: 6 species, North Sea: 9 species
cvA3yr <- matrix(NA,28,2)

lag=3

# Arrowtooth flounder
B <- arrowtoothFnum
idxF <- match(arrowtoothF$YEAR,exploitation$YEAR)
idxSST <- match(arrowtoothF$YEAR,sst$YEAR)

A <- data.frame(cbind(exploitation$Arrowtooth.flounder[idxF],sst$SST[idxSST-lag]))
M <- lifehist$M[1]

meanA3yr[1,] <- c(mean(A$X1/M),mean(A$X2))
cvA3yr[1,] <- c(sd(A$X1/M)/mean(A$X1/M),sd(A$X2)/mean(A$X2))

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
A <- data.frame(cbind(exploitation$Chilipepper.rockfish[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[2]

meanA3yr[2,] <- c(mean(A$X1/M),mean(A$X2))
cvA3yr[2,] <- c(sd(A$X1/lifehist$M[2])/mean(A$X1/lifehist$M[2]),sd(A$X2)/mean(A$X2))

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
A <- data.frame(cbind(exploitation$Dark.blotched.rockfish[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[3]

meanA3yr[3,] <- c(mean(A$X1/M),mean(A$X2))
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
A <- data.frame(cbind(exploitation$Dover.Sole[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[4]

meanA3yr[4,] <- c(mean(A$X1/M),mean(A$X2))
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
A <- data.frame(cbind(exploitation$English.Sole[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[5]

meanA3yr[5,] <- c(mean(A$X1/M),mean(A$X2))
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
A <- data.frame(cbind(exploitation$Lingcod[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[6]

meanA3yr[6,] <- c(mean(A$X1/M),mean(A$X2))
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
A <- data.frame(cbind(exploitation$Longspine.thornyhead[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[7]

meanA3yr[7,] <- c(mean(A$X1/M),mean(A$X2))
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
A <- data.frame(cbind(exploitation$Petrale.sole[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[8]

meanA3yr[8,] <- c(mean(A$X1/M),mean(A$X2))
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
idxF <- match(sardineS1$YEAR,exploitation$YEAR)
idxSST <- match(sardineS1$YEAR,sst$YEAR)
A <- data.frame(cbind(exploitation$Sardine[idxF],sst$SST[idxSST-lag]))[-4,]

M <- lifehist$M[9]

meanA3yr[9,] <- c(mean(A$X1/M),mean(A$X2))
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
A <- data.frame(cbind(exploitation$Splitnose.rockfish[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[10]

meanA3yr[10,] <- c(mean(A$X1/M),mean(A$X2))
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
A <- data.frame(cbind(exploitation$Yelloweye.rockfish[idxF],sst$SST[idxSST-lag]))

M <- lifehist$M[11]

meanA3yr[11,] <- c(mean(A$X1/M),mean(A$X2))
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

#AI pollock
meanA3yr[12,] <- NA
cvA3yr[12,] <- NA
adjR2Bout <- rbind(adjR2Bout,NA)
sigout <- rbind(sigout,NA)

#GOA pollock 
meanA3yr[13,] <- NA
cvA3yr[13,] <- NA
adjR2Bout <- rbind(adjR2Bout,NA)
sigout <- rbind(sigout,NA)

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

M <- lifehist$M[14]

meanA3yr[14,] <- c(mean(A$X1/M),mean(A$X2))
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
meanA3yr[15,] <- NA
cvA3yr[15,] <- NA
adjR2Bout <- rbind(adjR2Bout,NA)
sigout <- rbind(sigout,NA)

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

M <- lifehist$M[16]

meanA3yr[16,] <- c(mean(A$X1/M),mean(A$X2))
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

#GOA pacific cod
meanA3yr[17,] <- NA
cvA3yr[17,] <- NA
adjR2Bout <- rbind(adjR2Bout,NA)
sigout <- rbind(sigout,NA)

#GOA rex sole
meanA3yr[18,] <- NA
cvA3yr[18,] <- NA
adjR2Bout <- rbind(adjR2Bout,NA)
sigout <- rbind(sigout,NA)

# North Sea
year<-1977:2014
NSavg_tempQ1 <- read.csv("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/NSavg_tempQ1.csv",header = TRUE)

#cod
idxF <- match(year,cod$Year)
idxbottomT <- match(year,NSavg_tempQ1$YEAR)
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(cod$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(cod$Year[idxF],codLF$Year)
B <- codLF[idxB,-1] #LF w/o year

M <- cod$M[idxF]

meanA3yr[19,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[19,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(haddock$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(haddock$Year[idxF],haddockLF$Year)
B <- haddockLF[idxB,-1] #LF w/o year

M <- haddock$M[idxF]

meanA3yr[20,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[20,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(herring$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(herring$Year[idxF],herringLF$Year)
B <- herringLF[idxB,-1] #LF w/o year

M <- herring$M[idxF]

meanA3yr[21,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[21,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[!is.na(idxF)]
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[!is.na(idxbottomT)]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(mackerel$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(mackerel$Year[idxF],mackerelLF$Year)
B <- mackerelLF[idxB,-1] #LF w/o year

M <- lifehist$M[22]

meanA3yr[22,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[22,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[!is.na(idxF)]
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[!is.na(idxbottomT)]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(norwaypout$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(norwaypout$Year[idxF],norwaypoutLF$Year)
B <- norwaypoutLF[idxB,-1] #LF w/o year

M <- norwaypout$M[idxF]

meanA3yr[23,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[23,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(plaice$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(plaice$Year[idxF],plaiceLF$Year)
B <- plaiceLF[idxB,-1] #LF w/o year

M <- lifehist$M[24]

meanA3yr[24,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[24,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(saithe$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(saithe$Year[idxF],saitheLF$Year)
B <- saitheLF[idxB,-1] #LF w/o year

M <- lifehist$M[25]

meanA3yr[25,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[25,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(sole$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

idxB<- match(sole$Year[idxF],soleLF$Year)
B <- soleLF[idxB,-1] #LF w/o year

M <- lifehist$M[26]

meanA3yr[26,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[26,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(sprat$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

M <- sprat$M[idxF]

idxB<- match(sprat$Year[idxF],spratLF$Year)
B <- spratLF[idxB,-1] #LF w/o year

meanA3yr[27,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[27,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

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
idxF <- idxF[!is.na(idxF)]
idxF <- idxF[seq(1+lag,length(idxF))]
idxbottomT <- idxbottomT[!is.na(idxbottomT)]
idxbottomT <- idxbottomT[seq(1,length(idxbottomT)-lag)]

A <- data.frame(cbind(whiting$F[idxF],NSavg_tempQ1$temperature[idxbottomT]))
colnames(A)=c("exploitation","temperature")

M <- whiting$M[idxF[!is.na(idxF)]]

idxB<- match(whiting$Year[idxF],whitingLF$Year)
B <- whitingLF[idxB,-1] #LF w/o year

meanA3yr[28,] <- c(mean(A$exploitation/M),mean(A$temperature))
cvA3yr[28,] <- c(sd(A$exploitation/M)/mean(A$exploitation/M),sd(A$temperature)/mean(A$exploitation))

partB<-varpart(B, X= ~exploitation,~temperature,data=A)
plot(partB)
adjR2B<-c(partB$part$indfract$Adj.R.squared[1:3],partB$part$fract$Adj.R.squared[3])
B.rda <- rda(B~.,data=A)
plot(B.rda)
sig<-c(anova(rda(B~exploitation+Condition(temperature),data=A))$`Pr(>F)`[1],anova(rda(B~temperature+Condition(exploitation),data=A))$`Pr(>F)`[1])

adjR2Bout <- rbind(adjR2Bout,adjR2B)
sigout <- rbind(sigout,sig)

adjR2Bout[adjR2Bout<0]=0

write.csv(adjR2Bout,file = "output/adjR2Bout3yrTemp.csv")
write.csv(sigout,file = "output/sigout3yrTemp.csv")
write.csv(meanA3yr,file="output/meanSST3yr_F.csv")
write.csv(cvA3yr,file="output/cvSST3yr_F.csv")

rm(adjR2Bout,sigout)