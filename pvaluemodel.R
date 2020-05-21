find.threshold<-function(threshold,dataset){
  x = dataset[,1]
  y = dataset[,2]
  result = sum((x<threshold)*(y-0.5)^2)+sum((x>=threshold)*(y)^2) 
  return(result)
}

cal.pvalue<-function(tau,m,mean,sd){
  q = sqrt(m)*(mean-tau)/sd
  pvalue = (1-pnorm(q))
  return(pvalue)
}

find.tau<-function(tau,df,mean,sd){
  q = sqrt(df)*(mean-tau)/sd
  pvalue = (1-pnorm(q))
  result = sum((pvalue-0.5)^2)
  return(pvalue)
}

pvalue.model<-function(dataset){
  x=dataset[,1]
  y=dataset[,2]
  
  optionx=unique(x)
  m=mean(table(x))
  g.mean=sapply(1:length(optionx), function(i) mean(y[x==optionx[i]]))
  
  ybar.i = rep(g.mean,each=m)
  
  pooled.s = sqrt(sum((y-ybar.i)^2)/(length(y)-m))
  
  result=ga(type="real-valued", fitness=function(i) -find.tau(tau=i,df=m,mean=g.mean,sd=pooled.s), min=c(0), max=c(1))
  tau0=as.numeric(summary(result)$solution)
  
  pvalues=cal.pvalue(tau=tau0,m=m,mean=g.mean,sd=pooled.s)
  
  dataset2<-cbind(optionx,pvalues)
  
  result1=ga(type="real-valued",fitness=function(x) -find.threshold(threshold=x,dataset=dataset2), min=c(0), max=c(1))
  d0=as.numeric(summary(result1)$solution)[1]
  
  return(d0)
}
