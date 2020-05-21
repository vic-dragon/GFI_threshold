cal.x<-function(x){
  n=length(x)
  res=n*sum(x^2)-sum(x)^2
  return(res)
}

cal.xy<-function(x,y){
  n=length(x)
  res=n*sum(x*y)-sum(x)*sum(y)
  return(res)
}

cal.rss<-function(x,y){
  res=(cal.x(y)-cal.xy(x,y)^2/cal.x(x))/length(x)
  return(res)
}

f.dense.model.selection<-function(x, y, threshold, lambda){
  
  n=length(y)
  
  elambda<-exp(-lambda*(x-threshold))
  
  x1<-(1+elambda)^(-1)
  x2<-(x-threshold)*(1+elambda)^(-1)
  x3<-(x-threshold)^2*(1+elambda)^(-1)
  x4<-(x-threshold)^3*(1+elambda)^(-1)
  
  f.dense1= (-(n-1)/2)*log(cal.rss(x1,y)/2) 
  f.dense2= (-(n-1)/2)*log(cal.rss(x2,y)/2) 
  f.dense3= (-(n-1)/2)*log(cal.rss(x3,y)/2)
  f.dense4= (-(n-1)/2)*log(cal.rss(x4,y)/2)

  model.dense=c(f.dense1, f.dense2, f.dense3, f.dense4)
  
  co=round(max(model.dense))
    
  model.dense1=model.dense-co
    
  res=exp(model.dense1)/sum(exp(model.dense1))
  
  return(list(d1=res[1],d2=res[2],d3=res[3], d4=res[4]))
}