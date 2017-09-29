# Rearrange the temperature time series in Alaska
# 20161108 Crystal 
library(RNetCDF)
library(dplyr)
setwd("~/Desktop/fishsize_varpart/data/Alaska/")

#Aleutian Island annual SST
#Surface Gauss SST(C)
#Latitude Range used:   54.3 to  50.5
#Longitude Range used: 187.5 to 170.6
#NCEP Reanalysis 
#Produced at NOAA/ESRL PSD at http://www.esrl.noaa.gov/psd/data/timeseries/ 
#Climatology, if used, is now 1981-2010 (was 1971-2000). 
#Date submitted: 4/1/2015 at 02:12
AIsst <- read.table("AIannualSST.txt",header = TRUE)
colnames(AIsst) <- c("YEAR","SST")

#Gulf of Alaska 
GAK1 <- read.csv("GAK1.csv",header = TRUE)
GAK1surface <- filter(GAK1,Depth==30)
GAK1bottom <- filter(GAK1,Depth==150)
GAK1surface$YEAR <- round(GAK1surface$YEAR)
GAK1surface <- group_by(GAK1surface,YEAR)
GAK1surfaceT <- summarise(GAK1surface,T=mean(Temp),S=mean(Sal))
GAK1bottom$YEAR <- round(GAK1bottom$YEAR)
GAK1bottom <- group_by(GAK1bottom,YEAR)
GAK1bottomT <- summarise(GAK1bottom,T=mean(Temp),S=mean(Sal))

#EBS
nc<-open.nc("EBS/5820c0eea1d52-oc_Bottom.nc")
EBSbottom <- read.nc(nc,unpack = TRUE) #1982-2013
nc<-open.nc("EBS/5820c0eea1d52-oc_MaySST.nc") 
EBSsummerSST <- read.nc(nc,unpack = TRUE) #1948-2013
nc<-open.nc("EBS/5820c0eea1d52-oc_TsfcM2.nc")
EBSsurfactM2 <- read.nc(nc,unpack = TRUE) #1950-2013

YEAR=1948:2013
MaySST<-data.frame(cbind(YEAR,SST=EBSsummerSST$oc_MaySST_SST))
YEAR=1982:2013
summerBottomT <- data.frame(cbind(YEAR,summerBottomT=EBSbottom$oc_Bottom_BtmTmp))
YEAR=1950:2013
surfacTM2 <- data.frame(cbind(YEAR,sufaceT=EBSsurfactM2$oc_TsfcM2_TSFC))

#save(list = ls(all.names = TRUE),file="alaskaTemp.RData")
