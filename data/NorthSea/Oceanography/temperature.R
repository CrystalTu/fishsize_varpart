# North Sea near-bottom temperature and salinity (1977-)
# 20170121 Crystal

detach("package:plyr", unload=TRUE) #unload plyr and reshape2 to avoid conflict of "summarise" function
detach("package:reshape2",unload=TRUE)

library(dplyr)
library(tidyr)

setwd("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/")
rawdata <- read.csv("NS_TSdata.csv",header = TRUE)

# Extract month: Q1- 12,1,2 Q3-6,7,8
data <- select(rawdata,Cruise,Station,Year,Month,Bot..Depth..m.,TEMP..deg.C.)
temp <- filter(data,!is.na(Bot..Depth..m.)&!is.na(TEMP..deg.C.)) # exclude invalid data point
tempQ1 <- filter(temp,Month==12 | Month==1 | Month==2) 
tempQ3 <- filter(temp,Month==6 | Month==7 | Month==8) 

# Seasonal average for winter (Q1) and summer (Q3) 
tempQ1<-group_by(tempQ1,Year)
tempQ3<-group_by(tempQ3,Year)
NSavg_tempQ1<-summarise(tempQ1,depth=mean(Bot..Depth..m.),temperature=mean(TEMP..deg.C.))
NSavg_tempQ1<-mutate(NSavg_tempQ1,Quarter=1)
NSavg_tempQ3<-summarise(tempQ3,depth=mean(Bot..Depth..m.),temperature=mean(TEMP..deg.C.))
NSavg_tempQ3<-mutate(NSavg_tempQ3,Quarter=3)

# Combined into annual time series with mean and variation
#avg_temp <- bind_rows(avg_tempQ1,avg_tempQ3)
#avg_temp <- avg_temp[order(avg_temp$Year),]
#avg_temp <- group_by(avg_temp,Year)
#avg_temp <- summarise(avg_temp,meanDepth=mean(depth),meanTemp=mean(temperature),varTemp=var(temperature))

rm(rawdata,data,temp,tempQ1,tempQ3)
