# Re-arrange the IBTS length-structure data
# 20170117 Crystal
# 20170310 Crystal, full the absent length class with zero

#detach("package:plyr", unload=TRUE) #unload plyr and reshape2 to avoid conflict of "summarise" function
#detach("package:reshape2",unload=TRUE)

library(tidyr)
library(dplyr)

setwd("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/")
rawdata <- read.csv("CPUE per length per subarea_2016-09-21_1965to2016_q1q3.csv",header = TRUE)
targetSpecies <- c("Clupea harengus","Gadus morhua","Melanogrammus aeglefinus","Merlangius merlangus","Pleuronectes platessa","Pollachius virens","Scomber scombrus","Solea solea","Sprattus sprattus","Trisopterus esmarkii")

# trim year (1977-)
rawdata1977 <- filter(rawdata,Year>1976)

# seperate the target species
herring <- filter(rawdata1977,Species==targetSpecies[1])
cod <- filter(rawdata1977,Species==targetSpecies[2])
haddock <- filter(rawdata1977,Species==targetSpecies[3])
whiting <- filter(rawdata1977,Species==targetSpecies[4])
plaice <- filter(rawdata1977,Species==targetSpecies[5])
saithe <- filter(rawdata1977,Species==targetSpecies[6])
mackerel <- filter(rawdata1977,Species==targetSpecies[7])
sole <- filter(rawdata1977,Species==targetSpecies[8])
sprat <- filter(rawdata1977,Species==targetSpecies[9])
norwaypout <- filter(rawdata1977,Species==targetSpecies[10])

# winter total (Q1)  and length frequency
test<-filter(cod,Quarter==1)
test<-group_by(test,Year,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))

codLF<-spread(test,LngtClas,CPUE,fill=0) 

interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  codLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
codLF <- codLF[,c("Year",as.character(freq))] #sort the column with freq
codLF <- codLF[,-2] #remove 0mm

#haddock
test<-group_by(haddock,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
haddockLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  haddockLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
haddockLF <- haddockLF[,c("Year",as.character(freq))] #sort the column with freq
haddockLF <- haddockLF[,-2] #remove 0mm

#herring
test<-group_by(herring,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
herringLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 5
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  herringLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
herringLF <- herringLF[,c("Year",as.character(freq))] #sort the column with freq
herringLF <- herringLF[,-2] #remove 0mm

#whiting
test<-group_by(whiting,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
whitingLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  whitingLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
whitingLF <- whitingLF[,c("Year",as.character(freq))] #sort the column with freq
whitingLF <- whitingLF[,-2] #remove 0mm

#plaice
test<-group_by(plaice,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour)) #quarter sum
test<-filter(test,Quarter==1)
plaiceLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  plaiceLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
plaiceLF <- plaiceLF[,c("Year",as.character(freq))] #sort the column with freq
plaiceLF <- plaiceLF[,-2] #remove 0mm

#saithe
test<-group_by(saithe,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
saitheLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  saitheLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
saitheLF <- saitheLF[,c("Year",as.character(freq))] #sort the column with freq
saitheLF <- saitheLF[,-2] #remove 0mm

#mackerel
test<-group_by(mackerel,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
mackerelLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  mackerelLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
mackerelLF <- mackerelLF[,c("Year",as.character(freq))] #sort the column with freq
mackerelLF <- mackerelLF[,-2] #remove 0mm

#sole
test<-group_by(sole,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
soleLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  soleLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
soleLF <- soleLF[,c("Year",as.character(freq))] #sort the column with freq
soleLF <- soleLF[,-2] #remove 0mm

#sprat
test<-group_by(sprat,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
spratLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 5
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  spratLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
spratLF <- spratLF[,c("Year",as.character(freq))] #sort the column with freq
spratLF <- spratLF[,-2] #remove 0mm

#norwaypout
test<-group_by(norwaypout,Year,Quarter,LngtClas)
test<-summarise(test,CPUE=sum(CPUE_number_per_hour))
test<-filter(test,Quarter==1)
norwaypoutLF<-spread(test,LngtClas,CPUE,fill=0)
interval <- 10
Lmin <- min(unique(test$LngtClas))  
Lmax <- max(unique(test$LngtClas)) 
freq <- seq(Lmin,Lmax,by=interval)
absentClas<-match(freq,sort(unique(test$LngtClas))) #find the absent length class
xZero <- is.na(absentClas)
xZero <- as.character(freq[xZero])

for (n in 1:length(xZero)){
  norwaypoutLF[[xZero[n]]] <- 0 #fill-up the absent length class with zero
}
norwaypoutLF <- norwaypoutLF[,c("Year",as.character(freq))] #sort the column with freq
norwaypoutLF <- norwaypoutLF[,-2] #remove 0mm

rm(test,rawdata,rawdata1977)