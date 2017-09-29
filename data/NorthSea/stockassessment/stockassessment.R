# ICES Stock assessment table
# Note: annual fishing mortality from catch-weighted mean 
# of age-specific mortality

# 20170119 Crystal
# 20170817 Add natural mortality

library(XML)
library(plyr)
library(dplyr)
library(reshape2)
setwd("~/work/Fishery/ICES/stockassessment/")
source("naturalmortality/naturalmortality.R")

# 20170121 cod
test<-ldply(xmlToList("cod8052.xml"),data.frame)
cod <- filter(test,.id=="FishData")
cod<-test[-c(1:15),-c(1,2)]
idx <- sapply(cod,is.factor)
cod[,idx]<-lapply(cod[,idx],function(x) as.numeric(as.character(x)))
cod <- merge(cod,codM,by="Year")

# haddock
test<-ldply(xmlToList("haddock8068.xml"),data.frame)
haddock<-filter(test,.id=="FishData")
haddock <- haddock[,-c(1,2,20:30)]
idx <- sapply(haddock,is.factor)
haddock[,idx]<-lapply(haddock[,idx],function(x) as.numeric(as.character(x)))
haddock <- merge(haddock,haddockM,by="Year")

# herring
test<-ldply(xmlToList("herringAutumn7689.xml"),data.frame)
herring <- filter(test,.id=="FishData")
herring<-herring[,-c(1,2)]
idx <- sapply(herring,is.factor)
herring[,idx]<-lapply(herring[,idx],function(x) as.numeric(as.character(x)))
herring <- merge(herring,herringM,by="Year")

# mackerel
test<-ldply(xmlToList("mackerel8120.xml"),data.frame)
mackerel<-filter(test,.id=="FishData")
mackerel<-mackerel[,-c(1,2,19:42)]
idx <- sapply(mackerel,is.factor)
mackerel[,idx]<-lapply(mackerel[,idx],function(x) as.numeric(as.character(x)))

# norwaypout 
test<-ldply(xmlToList("norwaypout7998.xml"),data.frame)
norwaypout<-filter(test,.id=="FishData")
norwaypout <- norwaypout[-c(1,2)]
idx <- sapply(norwaypout,is.factor)
norwaypout[,idx]<-lapply(norwaypout[,idx],function(x) as.numeric(as.character(x)))
norwaypout <- merge(norwaypout,norwaypoutM,by="Year")

# plaice
test<-ldply(xmlToList("plaice7445.xml"),data.frame)
plaice<-filter(test,.id=="FishData")
plaice <- plaice[,-c(1,2,12:18)]
idx <- sapply(plaice,is.factor)
plaice[,idx]<-lapply(plaice[,idx],function(x) as.numeric(as.character(x)))

# saithe
test<-ldply(xmlToList("saithe8066.xml"),data.frame)
saithe<-filter(test,.id=="FishData")
saithe <- saithe[,-c(1,2,19:23)]
idx <- sapply(saithe,is.factor)
saithe[,idx]<-lapply(saithe[,idx],function(x) as.numeric(as.character(x)))

# sole
test<-ldply(xmlToList("sole7722.xml"),data.frame)
sole<-filter(test,.id=="FishData")
sole <- sole[,-c(1,2)]
idx <- sapply(sole,is.factor)
sole[,idx]<-lapply(sole[,idx],function(x) as.numeric(as.character(x)))

# sprat
test<-ldply(xmlToList("sprat7181.xml"),data.frame)
sprat <- filter(test,.id=="FishData")
sprat <- sprat[,-c(1,2)]
idx <- sapply(sprat,is.factor)
sprat[,idx]<-lapply(sprat[,idx],function(x) as.numeric(as.character(x)))
sprat <- merge(sprat,spratM,by="Year")

# whiting
test<-ldply(xmlToList("whiting7483.xml"),data.frame)
whiting<-filter(test,.id=="FishData")
whiting <- whiting[,-c(1,2,15:25)]
idx <- sapply(whiting,is.factor)
whiting[,idx]<-lapply(whiting[,idx],function(x) as.numeric(as.character(x)))
whiting <- merge(whiting,whitingM,by="Year")

rm(test,haddockM,herringM,norwaypoutM,spratM)
