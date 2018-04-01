# Rearrange the temperature time series in Alaska
# 20180212 Crystal 
library(dplyr)
setwd("~/Desktop/fishsize_varpart/data/Alaska/envi/")

#Aleutian triennial survey bottom and surface temperature
ai1983_2000 <- read.csv("ai1983_2000.csv")
ai2002_2012 <- read.csv("ai2002_2012.csv",header=TRUE,dec=",")

ai <- rbind(ai1983_2000,ai2002_2012,stringsAsFactors=F)
ai_pollock <- filter(ai, COMMON == "walleye pollock")
ai_pollock$BOT_TEMP[ai_pollock$BOT_TEMP==-9999.0] <- NA
ai_pollock$SURF_TEMP[ai_pollock$SURF_TEMP==-9999.0] <- NA
ai_pollock$SURF_TEMP[ai_pollock$NUMCPUE==-9999.0] <- NA
ai_pollock[,c(1,2,7,8,12,13,14)] <- sapply(ai_pollock[,c(1,2,7,8,12,13,14)],as.numeric)

ai_pollockAnnual <- group_by(ai_pollock,YEAR)
ai_pollockAnnual <- summarise(na.omit(ai_pollockAnnual),sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))
#ai1983_2000pollock <- group_by(ai1983_2000pollock,YEAR,STATION)
#ai1983_2000pollock <- summarise(ai1983_2000pollock,sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))

rm(ai1983_2000,ai2002_2012)

#EBS survey bottom and surface temperature
ebs1982_1984 <- read.csv("ebs1982_1984.csv")
ebs1985_1989 <- read.csv("ebs1985_1989.csv")
ebs1990_1994 <- read.csv("ebs1990_1994.csv")
ebs1995_1999 <- read.csv("ebs1995_1999.csv")
ebs2000_2004 <- read.csv("ebs2000_2004.csv")
ebs2005_2008 <- read.csv("ebs2005_2008.csv")
ebs2009_2012 <- read.csv("ebs2009_2012.csv")

ebs <- rbind(ebs1982_1984,ebs1985_1989,ebs1990_1994,ebs1995_1999,ebs2000_2004,ebs2005_2008,ebs2009_2012,stringsAsFactors=F)
ebs$BOT_TEMP[ebs$BOT_TEMP==-9999.0] <- NA
ebs$SURF_TEMP[ebs$SURF_TEMP==-9999.0] <- NA
ebs$SURF_TEMP[ebs$NUMCPUE==-9999.0] <- NA

ebs_pacificCod <- filter(ebs,COMMON == "Pacific cod")
ebs_pacificCod[,c(1,2,7,8,12,13,14)] <- sapply(ebs_pacificCod[,c(1,2,7,8,12,13,14)],as.numeric)
ebs_pacificCodAnnual <- group_by(ebs_pacificCod,YEAR)
ebs_pacificCodAnnual <- summarise(na.omit(ebs_pacificCodAnnual),sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))

ebs_flatheadSole <- filter(ebs,COMMON == "flathead sole")
ebs_flatheadSole[,c(1,2,7,8,12,13,14)] <- sapply(ebs_flatheadSole[,c(1,2,7,8,12,13,14)],as.numeric)
ebs_flatheadSoleAnnual <- group_by(ebs_flatheadSole,YEAR)
ebs_flatheadSoleAnnual <- summarise(na.omit(ebs_flatheadSoleAnnual),sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))

rm(ebs1982_1984,ebs1985_1989,ebs1990_1994,ebs1995_1999,ebs2000_2004,ebs2005_2008)

#GOA survey bottom and surface temperature
goa1984_1987 <- read.csv("goa1984_1987.csv")
goa1990_1999 <- read.csv("goa1990_1999.csv")
goa2001_2005 <- read.csv("goa2001_2005.csv")
goa2007_2013 <- read.csv("goa2007_2013.csv")

goa <- rbind(goa1984_1987,goa1990_1999,goa2001_2005,goa2007_2013)
goa$BOT_TEMP[goa$BOT_TEMP==-9999.0] <- NA
goa$SURF_TEMP[goa$SURF_TEMP==-9999.0] <- NA
goa$SURF_TEMP[goa$NUMCPUE==-9999.0] <- NA
goa[,c(1,2,7,8,12,13,14)] <- sapply(goa[,c(1,2,7,8,12,13,14)],as.numeric)

goa_pollock <- filter(goa, COMMON == "walleye pollock")
goa_flatheadSole <- filter(goa, COMMON == "flathead sole")
goa_pacificCod <- filter(goa, COMMON == "Pacific cod")
goa_rexSole <- filter(goa, COMMON == "rex sole")

goa_pollockAnnual <- group_by(goa_pollock,YEAR)
goa_pollockAnnual <- summarise(na.omit(goa_pollockAnnual),sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))
goa_flatheadSoleAnnual <- group_by(goa_flatheadSole,YEAR)
goa_flatheadSoleAnnual <- summarise(na.omit(goa_flatheadSoleAnnual),sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))
goa_pacificCodAnnual <- group_by(goa_pacificCod,YEAR)
goa_pacificCodAnnual <- summarise(na.omit(goa_pacificCodAnnual),sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))
goa_rexSoleAnnual <- group_by(goa_rexSole,YEAR)
goa_rexSoleAnnual <- summarise(na.omit(goa_rexSoleAnnual),sum(WTCPUE),sum(NUMCPUE),mean(BOT_DEPTH),mean(BOT_TEMP),mean(SURF_TEMP))

rm(goa1984_1987,goa1990_1999,goa2001_2005,goa2007_2013)

colnames(ai_pollockAnnual)<-c("YEAR","WTCPUE","NUMCPUE","BOT_DEPTH","BOT_TEMP","SURF_TEMP")
colnames(ebs_pacificCodAnnual)<-c("YEAR","WTCPUE","NUMCPUE","BOT_DEPTH","BOT_TEMP","SURF_TEMP") 
colnames(ebs_flatheadSoleAnnual)<-c("YEAR","WTCPUE","NUMCPUE","BOT_DEPTH","BOT_TEMP","SURF_TEMP") 
colnames(goa_flatheadSoleAnnual)<-c("YEAR","WTCPUE","NUMCPUE","BOT_DEPTH","BOT_TEMP","SURF_TEMP") 
colnames(goa_pacificCodAnnual)<-c("YEAR","WTCPUE","NUMCPUE","BOT_DEPTH","BOT_TEMP","SURF_TEMP") 
colnames(goa_pollockAnnual)<-c("YEAR","WTCPUE","NUMCPUE","BOT_DEPTH","BOT_TEMP","SURF_TEMP") 
colnames(goa_rexSoleAnnual)<-c("YEAR","WTCPUE","NUMCPUE","BOT_DEPTH","BOT_TEMP","SURF_TEMP") 

write.csv(ai_pollockAnnual,file="ai_pollockAnnual.csv",row.names = FALSE)
write.csv(ebs_pacificCodAnnual,file="ebs_pacificCodAnnual.csv",row.names = FALSE)
write.csv(ebs_flatheadSoleAnnual,file="ebs_flatheadSoleAnnual.csv",row.names = FALSE)
write.csv(goa_flatheadSoleAnnual,file="goa_flatheadSoleAnnual.csv",row.names = FALSE)
write.csv(goa_pacificCodAnnual,file="goa_pacificCodAnnual.csv",row.names = FALSE)
write.csv(goa_pollockAnnual,file="goa_pollockAnnual.csv",row.names = FALSE)
write.csv(goa_rexSoleAnnual,file="goa_rexSoleAnnual.csv",row.names = FALSE)

rm(goa,ai,ebs)
