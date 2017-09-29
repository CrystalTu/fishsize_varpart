# Indicator function
# Crystal Tu 2016/12/06

fishsize <- function(x,sizebin){
  result <- matrix(rep(sizebin,each=nrow(x)),nrow(x),ncol(x))
  return(result)
}

L95perc <- function(x,fishsize){
  out<-rep(0,nrow(x))
  perc<-apply(x,1,quantile,probs=.95)
  for (i in 1:nrow(x)){
    idx<-which((abs(x[i,]-perc[i]))==min(abs(x[i,]-perc[i])))
#    idx<-c(idx-2,idx-1,idx,idx+1,idx+2)
#    out[i]<-interp1(as.numeric(x[i,idx]),sizebin[idx],perc[i],method="spline")
  out[i]<-as.numeric(max(fishsize[i,idx]))  
  }
  result<-out
}

meansize <- function(x,fishsize){
  result <- rowSums(fishsize*x)/rowSums(x)
  return(result)
}

shannon <- function(x){
  result <- rowSums(x*log(x),na.rm=TRUE)
  return(result)
}

evenness <- function(x){
  result <- shannon(x)/log(rowSums(ifelse(x>0,1,0)))
  return(result)
}

skewness <- function(x){
  mean_matrix <- matrix(rep(rowMeans(x),time=ncol(x)),nrow(x),ncol(x))
  std <- matrix(0,nrow(x),1)
  for (i in 1:nrow(x)) {
    std[i] <- sd(x[i,])
  }
  result <- (1/ncol(x))*rowSums(((x-mean_matrix)/std)^3)
  return(result)
}