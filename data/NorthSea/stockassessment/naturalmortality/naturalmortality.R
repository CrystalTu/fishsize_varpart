# Calculate the natural mortality (M) corresponding to the fishing mortality
# North Sea stock 
#
# Crystal Tu 2017/08/15

#cod
cod <- read.csv('naturalmortality/NScod.csv',header = TRUE)
codM <- data.frame(Year=cod$Year,M=rowMeans(cod[,3:5])) #F2-4

#haddock
haddock <- read.csv('naturalmortality/NShaddock.csv',header = TRUE)
haddockM <- data.frame(Year=haddock$Year,M=rowMeans(haddock[,3:5])) #F2-4

#herring
herring <- read.csv('naturalmortality/NSherring.csv',header = TRUE)
herringM <- data.frame(Year=herring$Year,M=rowMeans(herring[,4:8])) #F2-6

#norway pout
norwaypout <- read.csv('naturalmortality/NSnorwaypout.csv',header = TRUE)
norwaypoutM <- data.frame(Year=norwaypout$Year,M=rowMeans(norwaypout[,3:4])) #F1-2

#sprat 
sprat <- read.csv('naturalmortality/NSsprat.csv',header=TRUE)
spratM <- data.frame(Year=sprat$Year,M=rowMeans(sprat[,2:3])) #F1-2

#whiting
whiting <- read.csv('naturalmortality/NSwhiting.csv',header = TRUE)
whitingM <- data.frame(Year=whiting$Year,M=rowMeans(whiting[,3:7])) #F2-6