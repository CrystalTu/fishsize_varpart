# Pre-set of length distribution data for further analysis
# 20161004 Crystal
# 20180211 Use survey length composition in Alaska
library(dplyr)
library(tidyr)
detach("package:plyr", unload=TRUE)

setwd("~/Desktop/fishsize_varpart/data/WestUS/length/")

# West US
arrowtooth <- read.csv("ArrowtoothFlounder.csv",header = TRUE)
arrowtoothF <- filter(arrowtooth,gender==2)
arrowtoothFnum <- select(arrowtoothF,-c(YEAR,gender,Nsamp))*arrowtoothF$Nsamp

chilipepperRF <- read.csv("ChilipepperRockfish.csv",header=TRUE)
chilipepperRFf <- filter(chilipepperRF,gender==2)
chilipepperRFfnum <- select(chilipepperRFf,-c(YEAR,gender,Nsamp))*chilipepperRFf$Nsamp

darkblotchedRF <- read.csv("DarkBlotchRockfish.csv",header = TRUE)
darkblotchedRFf <- filter(darkblotchedRF,gender==2)
darkblotchedRFfnum <- select(darkblotchedRFf,-c(YEAR,fleet,seas,gender,part,gender,effn))*darkblotchedRFf$effn

dover <- read.csv("DoverSole.csv",header = TRUE) #use both sex
doverN <- filter(dover,region==1,gender==2) + filter(dover,region==1,gender==1)
doverN$YEAR <- 1966:2004
doverS <- filter(dover,region==2,gender==2,YEAR>1967) + filter(dover,region==2,gender==1,YEAR>1967)
doverS$YEAR <- 1969:2004
dovernum <- (doverN[-c(1:3),4:23]+doverS[,4:23])/2 #Dover total year == doverS

english <- read.csv("EnglishSole.csv",header=TRUE)
englishF <- filter(english,gender==2)
englishFnum <- select(englishF,-c(YEAR,seas,type,gender,partition,Nsamp))*englishF$Nsamp

lingcod <- read.csv("lingcod.csv",header=TRUE)
lingcodF <- filter(lingcod,gender==2)
lingcodFnum <- select(lingcodF,-c(YEAR,gender,part,effn))

longspine <- read.csv("LongspinThornyhead.csv",header=TRUE)
longspineNum <- select(longspine,-c(YEAR,agency,X.trip,X.fish))*longspine$X.fish

petraleSole <- read.csv("petraleSole.csv",header=TRUE)
petraleSoleFnum <- (petraleSole[,7:32]/100)*petraleSole$nSamps 

sardine <- read.csv("sardine.csv",header=TRUE)
sardineS1 <- filter(sardine,seas==1)
sardineS1num <- select(sardineS1,-c(YEAR,seas,flt.svy,gender,part,Nsamp))*sardineS1$Nsamp
sardineS2 <- filter(sardine,seas==2)
sardineS2num <- select(sardineS2,-c(YEAR,seas,flt.svy,gender,part,Nsamp))*sardineS2$Nsamp
sardineNum <- (sardineS1num + sardineS2num)/2 #missing data in 1984

splitnoseRF <- read.csv("SplitnoseRockfish.csv",header = TRUE)
splitnoseRFf <- filter(splitnoseRF,gender==2)
splitnoseRFfnum <- (select(splitnoseRFf,-c(YEAR,sea,flt.svy,gender,part,Nsamp))/100)*splitnoseRFf$Nsamp

yelloweyeRF <- read.csv("YelloweyeRockfish.csv",header = TRUE)
yelloweyeRFnum <- select(yelloweyeRF,-c(YEAR,seas,type,gender,partition,Nsamp))

#Alaska
setwd("~/Desktop/fishsize_varpart/data/Alaska/length/")

#Aleutian Island pollock into length-bin frequency
AIpollock <- read.csv("AIpollock.csv",header = TRUE)
annual <- group_by(AIpollock,YEAR,LENGTH_cm)
AIpollockLF <- spread(summarise(annual,total=sum(FREQUENCY)),LENGTH_cm,total)
AIpollockLF[is.na(AIpollockLF)] <- 0

#Insert zero vector in zero length bin
lengthMin <- min(AIpollock$LENGTH_cm)
lengthMax <- max(AIpollock$LENGTH_cm)
lengthBin <- lengthMin:lengthMax
bin <- as.numeric(unique(colnames(AIpollockLF))[-1])
instCol <- as.character(lengthBin[-match(bin,lengthBin)])
tmp <- data.frame(matrix(0,length(AIpollockLF$YEAR), length(instCol)[1]),stringsAsFactors = F)
tmp <- cbind(AIpollockLF$YEAR,tmp)
colnames(tmp) <- c("YEAR",instCol)
AIpollocknum <- left_join(AIpollockLF,tmp,by="YEAR")
AIpollocknum <- AIpollocknum[,match(c("YEAR",as.character(lengthBin)),names(AIpollocknum))]
AIpollocknum <- AIpollocknum[,-1]

#GOA pollock
GOApollock <- read.csv("GOApollock.csv",header=TRUE)
GOApollockNum <- GOApollock[,-1]

#BSAI Flathead sole
BSAIflatheadSole <- read.csv("BSAIflathead_sole.csv",header=TRUE)
BSAIflatheadSole[,c(3:26)] <- sapply(BSAIflatheadSole[,c(3:26)],as.numeric)
BSAIflatheadSoleF <- filter(BSAIflatheadSole,gender==2)[-1,]
BSAIflatheadSoleM <- filter(BSAIflatheadSole,gender==1)[-1,]
BSAIflatheadSoleNum <- (BSAIflatheadSoleF + BSAIflatheadSoleM)[,-c(1:2)]

#GOA Flathead sole
GOAflatheadSole <- read.csv("GOAflathead_sole.csv",header = TRUE)
GOAflatheadSole[,c(3:16)]<-sapply(GOAflatheadSole[,c(3:16)],as.numeric)
GOAflatheadSoleF <- filter(GOAflatheadSole,gender==2)
GOAflatheadSoleM <- filter(GOAflatheadSole,gender==1)
GOAflatheadSoleNum <- (GOAflatheadSoleF + GOAflatheadSoleM)[,-c(1:2)]

#EBS Pacific cod
EBSpacifcCod <- read.csv("EBSpacific_cod.csv",header = TRUE)
EBSpacifcCodNum <- EBSpacifcCod[,-1]

#GOA Pacific cod
GOApacificCod <- read.csv("GOApacific_cod.csv",header = TRUE)
GOApacificCodNum <- GOApacificCod[,-1]

#GOA Rex sole
GOArexSole <- read.csv("GOArexSole.csv",header = TRUE)
GOArexSole[,c(3:18)] <- sapply(GOArexSole[,c(3:18)],as.numeric) 
GOArexSoleF <- filter(GOArexSole,gender==2)
GOArexSoleM <- filter(GOArexSole,gender==1)
GOArexSoleNum <- (GOArexSoleF + GOArexSoleM)[,-c(1:2)]

rm(bin,instCol,lengthBin,lengthMax,lengthMin,tmp,annual)
