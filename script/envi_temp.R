# Temperature plot
# Crystal 20170306
setwd("~/Desktop/fishsize_varpart/")
setEPS()
postscript("output/FigS2_envi_temp.eps")
par(mfrow=c(3,2))

#West US
west_us_SST <- read.table("~/Desktop/fishsize_varpart/data/WestUS/WestUSannualSSTfrom1948.txt",header = FALSE)
colnames(west_us_SST)<-c("YEAR","SST")
YEAR<-1965:2008
idxT<-match(west_us_SST$YEAR,YEAR)
plot(west_us_SST$YEAR[idxT],west_us_SST$SST[idxT],type="l",main="West US", xlab="",ylab="")

#Alaska (dot plot for AI and GOA)
AItemp <- read.csv("data/Alaska/envi/ai_pollockAnnual.csv",header=TRUE)
plot(AItemp$YEAR,AItemp$BOT_TEMP,type="b",main="AI", xlab="",ylab="")
EBStemp <- read.csv("data/Alaska/envi/ebs_flatheadSoleAnnual.csv",header = TRUE)
plot(EBStemp$YEAR,EBStemp$BOT_TEMP,type="l",main="EBS", xlab="",ylab="")
GOAtemp <- read.csv("data/Alaska/envi/goa_pollockAnnual.csv",header=TRUE)
plot(GOAtemp$YEAR,GOAtemp$BOT_TEMP,type="b",main="GOA", xlab="",ylab="")
#North Sea
NSavg_tempQ1 <- read.csv(file = "~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/NSavg_tempQ1.csv",header = TRUE)
NSavg_tempQ1_surf <- read.csv(file = "~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/NSavg_tempQ1_surf.csv",header = TRUE)
plot(NSavg_tempQ1$YEAR,NSavg_tempQ1$temperature,type="l",main="North Sea", xlab="",ylab="")
dev.off()

#Surface water temperature
setEPS()
postscript("output/FigS3_envi_temp_surf.eps")
par(mfrow=c(2,2))

plot(AItemp$YEAR,AItemp$SURF_TEMP,type="b",main="AI", xlab="",ylab="")
plot(EBStemp$YEAR,EBStemp$SURF_TEMP,type="l",main="EBS", xlab="",ylab="")
plot(GOAtemp$YEAR,GOAtemp$SURF_TEMP,type="b",main="GOA", xlab="",ylab="")
plot(NSavg_tempQ1_surf$YEAR,NSavg_tempQ1_surf$temperature,type="l",main="North Sea", xlab="",ylab="")
dev.off()

#correlation
cor(cbind(AItemp$SURF_TEMP,AItemp$BOT_TEMP))
cor(cbind(EBStemp$SURF_TEMP,EBStemp$BOT_TEMP))
cor(cbind(GOAtemp$SURF_TEMP,GOAtemp$BOT_TEMP))
cor(cbind(NSavg_tempQ1_surf$temperature,NSavg_tempQ1$temperature))