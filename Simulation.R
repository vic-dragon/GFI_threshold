parameter.quantile<-function(result,lb=0.025,ub=0.975){
  result=rbind(quantile(result[,1],probs=c(lb,0.5,ub)),quantile(result[,2],probs=c(lb,0.5,ub)),
               quantile(result[,3],probs=c(lb,0.5,ub)),quantile(result[,4],probs=c(lb,0.5,ub)))
  rownames(result)<-c("alpha0","alpha1","sigma2","threshold")
  return(result)
}

estimated.parameter<-function(result){
  result=cbind(mean(result[,1]),mean(result[,2]),mean(result[,3]),mean(result[,4]))
  colnames(result)<-c("alpha0","alpha1","sigma2","threshold")
  return(result)
}


Fid.spline<-function(dataset, old.threshold, sd.threshold, burn.in=1000, s.size=2000){
  
  x=dataset[,1]
  y=dataset[,2]
  
  lambda=log(999)/(min(diff(unique(x))))*2
  
  case=Metro.alg.spline(x, y, old.threshold=old.threshold, sd.threshold = sd.threshold, burn.in=burn.in, s.size=s.size,lambda=lambda)
  loc.sample=1:(s.size/10)*10
  case=case[loc.sample,]
  
  est.para=estimated.parameter(case)
  CI.95=parameter.quantile(case)
  
  pvaluemethod=pvalue.model(dataset)
  
  psi.value=runif(1,min=min(x),max=max(x))
  lm_seg=lm(y~x)
  result.seg=tryCatch(segmented(lm_seg,seg.Z=~x,psi=psi.value), error=function(e) NA)
  #indi.seg=tryCatch(result.seg$seed, error=function(e) NA)
  
  #if(length(is.na(indi.seg))==0){
  if(length(result.seg)<3){
    segmethod=c(0,0,0)
  }else{
    seg.est=result.seg$psi[2]
    seg.ub=result.seg$psi[2]+1.96*result.seg$psi[3]
    seg.lb=result.seg$psi[2]-1.96*result.seg$psi[3]
    segmethod=c(seg.est,seg.lb,seg.ub)
  }
  names(segmethod)<-c("seg.est","seg.lb","seg.ub")

  return(list(model=table(case[,5]), est.para=est.para, CI=CI.95, pvaluemethod=pvaluemethod, segmethod=segmethod))
}

