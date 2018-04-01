# North Sea near-bottom temperature and salinity from ICES Oceanography dataset (1977-)
# 20170121 Crystal

#detach("package:plyr", unload=TRUE) #unload plyr and reshape2 to avoid conflict of "summarise" function
#detach("package:reshape2",unload=TRUE)

library(dplyr)
library(tidyr)
detach(package:plyr)

library(dplyr)
library(tidyr)
detach(package:plyr)

setwd("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/")
rawdata <- read.csv("NS_TSdata.csv",header = TRUE)

# Extract month: Q1- 1,2 Q3-7,8,9
data <- select(rawdata,Cruise,Station,Year,Month,Bot..Depth..m.,TEMP..deg.C.)
temp <- filter(data,!is.na(Bot..Depth..m.)&!is.na(TEMP..deg.C.)) # exclude invalid data point
tempQ1 <- filter(temp,Month==1 | Month==2) 
tempQ3 <- filter(temp,Month==7 | Month==8 | Month==9) 

# Seasonal average for winter (Q1) and summer (Q3) 
tempQ1<-group_by(tempQ1,Year)
tempQ3<-group_by(tempQ3,Year)
NSavg_tempQ1<-summarise(tempQ1,depth=mean(Bot..Depth..m.),temperature=mean(TEMP..deg.C.))
NSavg_tempQ1<-mutate(NSavg_tempQ1,Quarter=1)
NSavg_tempQ3<-summarise(tempQ3,depth=mean(Bot..Depth..m.),temperature=mean(TEMP..deg.C.))
NSavg_tempQ3<-mutate(NSavg_tempQ3,Quarter=3)

rm(rawdata,data,temp,tempQ1,tempQ3)
write.csv(NSavg_tempQ1,file = "NSavg_tempQ1.csv",row.names = FALSE)
write.csv(NSavg_tempQ3,file = "NSavg_tempQ3.csv",row.names = FALSE)