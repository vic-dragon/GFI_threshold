organize.result<-function(result){
  
  n=length(result)

  org.res=matrix(ncol=24,nrow=n)
  
  colnames(org.res)=c("model1","model2","model3","model4","mean.alpha0","alpha0lb","alpha0ub","mean.alpha1","alpha1lb","alpha1ub",
                      "mean.sigma","sigmalb","sigmaub","mean.threshold","thresholdlb","thresholdub","method2","med.alpha0","med.alpha1","med.sigma","med.threshold","seg.est","seg.lb","seg.ub")
  
  for(i in 1:n){
    
    c1=result[[i]]$model[names(result[[i]]$model)==1]
    if(length(c1)==0){c1=0}
    c2=result[[i]]$model[names(result[[i]]$model)==2]
    if(length(c2)==0){c2=0}
    c3=result[[i]]$model[names(result[[i]]$model)==3]
    if(length(c3)==0){c3=0}
    c4=result[[i]]$model[names(result[[i]]$model)==4]
    if(length(c4)==0){c4=0}
    
    org.res[i,1]<-c1
    org.res[i,2]<-c2
    org.res[i,3]<-c3
    org.res[i,4]<-c4
    
    org.res[i,5]<-result[[i]]$est.para[1]
    org.res[i,8]<-result[[i]]$est.para[2]
    org.res[i,11]<-result[[i]]$est.para[3]
    org.res[i,14]<-result[[i]]$est.para[4]
    
    org.res[i,6]<-result[[i]]$CI[1]
    org.res[i,9]<-result[[i]]$CI[2]
    org.res[i,12]<-result[[i]]$CI[3]
    org.res[i,15]<-result[[i]]$CI[4]
    
    org.res[i,7]<-result[[i]]$CI[9]
    org.res[i,10]<-result[[i]]$CI[10]
    org.res[i,13]<-result[[i]]$CI[11]
    org.res[i,16]<-result[[i]]$CI[12]
    
    org.res[i,17]<-result[[i]]$pvaluemethod
    
    org.res[i,18]<-result[[i]]$CI[5]
    org.res[i,19]<-result[[i]]$CI[6]
    org.res[i,20]<-result[[i]]$CI[7]
    org.res[i,21]<-result[[i]]$CI[8]
    
    org.res[i,22]<-result[[i]]$segmethod[1]
    org.res[i,23]<-result[[i]]$segmethod[2]
    org.res[i,24]<-result[[i]]$segmethod[3]
    
  }
  return(as.data.frame(org.res))
}

sum.result<-function(result,true.threshold=0.5){

  res=organize.result(result)
  n=dim(res)[1]
  selected.model=table(sapply(1:n, function(x) which.max(res[x,1:3])))
  
  fid.mean=mean(res$mean.threshold)
  fid.se=sd(res$mean.threshold)/sqrt(n)
  
  pval.mean=mean(res$method2)
  pval.se=sd(res$method2)/sqrt(n)
  
  mean.rmse=sqrt(sum((res$mean.threshold - true.threshold)^2)/n)
  med.rmse=sqrt(sum((res$med.threshold - true.threshold)^2)/n)
  pvalue.rmse=sqrt(sum((res$method2 - true.threshold)^2)/n)
  
  #test.res<-t.test(abs(res$mean.threshold-true.threshold)-abs(res$method2-true.threshold))
  #test.res$p.value
  #test.res$estimate
  
  range=mean(res$thresholdub-res$thresholdlb)
  cov.rate=sum(res$thresholdlb<true.threshold & res$thresholdub>true.threshold)/n
  
  seg.num<-sum(res$seg.est!=0)
  seg.rmse<-sqrt(sum((res$seg.est[res$seg.est!=0] - true.threshold)^2)/seg.num)
  seg.range=mean(res$seg.ub[res$seg.est!=0]-res$seg.lb[res$seg.est!=0])
  seg.rate=sum(res$seg.lb[res$seg.est!=0]<true.threshold & res$seg.ub[res$seg.est!=0]>true.threshold)/seg.num
  
  seg.mean=mean(res$seg.est[res$seg.est!=0])
  seg.se=sd(res$seg.est[res$seg.est!=0])/sqrt(seg.num)
  
  
  
  return(list(model=selected.model,  
              fid.mean=fid.mean, fid.se=fid.se, pval.mean=pval.mean, pval.se=pval.se,
              seg.mean=seg.mean, seg.se=seg.se,
              fid.rmse=round(mean.rmse,3), seg.rmse=round(seg.rmse,3),
              pvalue.rmse=round(pvalue.rmse,3), 
              cov.rate=round(cov.rate,3), range=round(range,3),
              seg.rate=round(seg.rate,3), seg.range=round(seg.range,3), seg.num=seg.num
              ))
  
  #return(list(model=selected.model,  
  #            mean.rmse=round(mean.rmse,3), med.rmse=round(med.rmse,3), 
  #            pvalue.rmse=round(pvalue.rmse,3), seg.rmse=round(seg.rmse,3),
  #            test.est=round(test.res$estimate,3), test.pvalue=round(test.res$p.value,3), 
  #            cov.rate=round(cov.rate,3), range=round(range,3),
  #            seg.rate=round(seg.rate,3), seg.range=round(seg.range,3),
  #            mean1=round(mean(res$mean.threshold),3), mean2=round(mean(res$method2),3)
  #            ))
}
to_vector<-function(result){
  return(c(result$fid.mean, result$fid.se, 
           result$pval.mean, result$pval.se,
           result$seg.mean, result$seg.se,
           result$fid.rmse, result$pvalue.rmse,
           result$seg.rmse, 
           result$cov.rate, result$range,
           result$seg.rate, result$seg.range,
           result$seg.num))
}

setwd('/Users/Seungyong/Desktop/Fiducial Threshold project/result')

#Model selection result for fiducial inference.
org_result_model=as.data.frame(matrix(0,ncol=25,nrow=11))
colnames(org_result_model)=c('setting',
                             'f_model1_f1','f_model2_f1','f_model3_f1','f_model4_f1',
                             'p_model1_f1','p_model2_f1','p_model3_f1','p_model4_f1',
                             'f_model1_f2','f_model2_f2','f_model3_f2','f_model4_f2',
                             'p_model1_f2','p_model2_f2','p_model3_f2','p_model4_f2',
                             'f_model1_f3','f_model2_f3','f_model3_f3','f_model4_f3',
                             'p_model1_f3','p_model2_f3','p_model3_f3','p_model4_f3')
n_values=c(6,6,10,10,10,10,20,20,50,50)
m_values=c(5,5,10,10,20,20,50,50,100,100)
s_values=c(1,3,1,3,1,3,1,3,1,3)

for (j in 1:10){
  file_name=paste0('n',n_values[j],'m',m_values[j],'s0',s_values[j])
  eval(parse(text = paste0('load("',file_name,'.Rdata")')))
  org_result_model[j,1] = file_name
  
  org_res1 = organize.result(store.f1)
  org_res2 = organize.result(store.f2)
  org_res3 = organize.result(store.f3)

  table_model=table(apply(org_res1[,1:4],1,which.max))/200
  org_result_model[j,2]=ifelse(length(is.na(which(names(table_model)=="1")))==0, 0, table_model[which(names(table_model)=="1")])
  org_result_model[j,3]=ifelse(length(is.na(which(names(table_model)=="2")))==0, 0, table_model[which(names(table_model)=="2")])
  org_result_model[j,4]=ifelse(length(is.na(which(names(table_model)=="3")))==0, 0, table_model[which(names(table_model)=="3")])
  org_result_model[j,5]=ifelse(length(is.na(which(names(table_model)=="4")))==0, 0, table_model[which(names(table_model)=="4")])
  
  org_result_model[j,6]=sum(org_res1$model1)/(40000)
  org_result_model[j,7]=sum(org_res1$model2)/(40000)
  org_result_model[j,8]=sum(org_res1$model3)/(40000)
  org_result_model[j,9]=sum(org_res1$model4)/(40000)
  
  table_model=table(apply(org_res2[,1:4],1,which.max))/200
  org_result_model[j,10]=ifelse(length(is.na(which(names(table_model)=="1")))==0, 0, table_model[which(names(table_model)=="1")])
  org_result_model[j,11]=ifelse(length(is.na(which(names(table_model)=="2")))==0, 0, table_model[which(names(table_model)=="2")])
  org_result_model[j,12]=ifelse(length(is.na(which(names(table_model)=="3")))==0, 0, table_model[which(names(table_model)=="3")])
  org_result_model[j,13]=ifelse(length(is.na(which(names(table_model)=="4")))==0, 0, table_model[which(names(table_model)=="4")])
  
  org_result_model[j,14]=sum(org_res2$model1)/(40000)
  org_result_model[j,15]=sum(org_res2$model2)/(40000)
  org_result_model[j,16]=sum(org_res2$model3)/(40000)
  org_result_model[j,17]=sum(org_res2$model4)/(40000)
  
  table_model=table(apply(org_res3[,1:4],1,which.max))/200
  org_result_model[j,18]=ifelse(length(is.na(which(names(table_model)=="1")))==0, 0, table_model[which(names(table_model)=="1")])
  org_result_model[j,19]=ifelse(length(is.na(which(names(table_model)=="2")))==0, 0, table_model[which(names(table_model)=="2")])
  org_result_model[j,20]=ifelse(length(is.na(which(names(table_model)=="3")))==0, 0, table_model[which(names(table_model)=="3")])
  org_result_model[j,21]=ifelse(length(is.na(which(names(table_model)=="4")))==0, 0, table_model[which(names(table_model)=="4")])
  
  org_result_model[j,22]=sum(org_res3$model1)/(40000)
  org_result_model[j,23]=sum(org_res3$model2)/(40000)
  org_result_model[j,24]=sum(org_res3$model3)/(40000)
  org_result_model[j,25]=sum(org_res3$model4)/(40000)
}


s01_model=org_result_model[c(1:5)*2-1,]
s03_model=org_result_model[c(1:5)*2,]

# 3 methods comparison.

org_result=as.data.frame(matrix(0,ncol=18,nrow=1))
colnames(org_result)=c('n','m','sigma','fnct','fid_mean','fid_se','pval_mean','pval_se','seg_mean','seg_se',
                       'fid_rmse','pval_rmse','seg_rmse','fid_cr','fid_range','seg_cr','seg_range','seg_num')

n_values=c(6,6,6,6,10,10,10,10,10,10,20,20,20,20,50,50)
m_values=c(5,5,10,10,5,5,10,10,20,20,10,10,50,50,100,100)
s_values=c(1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3)

org_result_pvalue=as.data.frame(matrix(0,ncol=25,nrow=11))
n_values=c(6,6,10,10,10,10,20,20,50,50)
m_values=c(5,5,10,10,20,20,50,50,100,100)
s_values=c(1,3,1,3,1,3,1,3,1,3)
names(org_result_pvalue)=c('setting','pval12_f1','pval13_f1','pval12_f2','pval13_f2','pval12_f3','pval13_f3',
                          'mse_gfi_f1','mse_mle_f1','mse_p_f1',
                          'mse_gfi_f2','mse_mle_f2','mse_p_f2',
                          'mse_gfi_f3','mse_mle_f3','mse_p_f3',
                          'se_gfi_f1','se_mle_f1','se_p_f1',
                          'se_gfi_f2','se_mle_f2','se_p_f2',
                          'se_gfi_f3','se_mle_f3','se_p_f3')
for (j in 1:10){
  file_name=paste0('n',n_values[j],'m',m_values[j],'s0',s_values[j])
  eval(parse(text = paste0('load("',file_name,'.Rdata")')))
  org_result_pvalue[j,1] = file_name

  org_res1 = organize.result(store.f1)
  m11=(org_res1$mean.threshold - 0.5)^2
  m12=(org_res1$seg.est - 0.5)^2
  m13=(org_res1$method2 - 0.5)^2
  
  test1_12=t.test(m11,m12,alternative='less',paired=T)
  test1_13=t.test(m11,m13,alternative='less',paired=T)
  
  org_result_pvalue[j,2] = test1_12$p.value
  org_result_pvalue[j,3] = test1_13$p.value
  
  org_res2 = organize.result(store.f2)
  m21=(org_res2$mean.threshold - 0.5)^2
  m22=(org_res2$seg.est - 0.5)^2
  m23=(org_res2$method2 - 0.5)^2
  test2_12=t.test(m21,m22,alternative='less',paired=T)
  test2_13=t.test(m21,m23,alternative='less',paired=T)
  org_result_pvalue[j,4] = test2_12$p.value
  org_result_pvalue[j,5] = test2_13$p.value
  
  org_res3 = organize.result(store.f3)
  m31=(org_res3$mean.threshold - 0.5)^2
  m32=(org_res3$seg.est - 0.5)^2
  m33=(org_res3$method2 - 0.5)^2
  test3_12=t.test(m31,m32,alternative='less',paired=T)
  test3_13=t.test(m31,m33,alternative='less',paired=T)
  org_result_pvalue[j,6] = test3_12$p.value
  org_result_pvalue[j,7] = test3_13$p.value
  
  org_result_pvalue[j,8] = sum(m11)/200
  org_result_pvalue[j,9] = sum(m12)/200
  org_result_pvalue[j,10] = sum(m13)/200
  
  org_result_pvalue[j,11] = sum(m21)/200
  org_result_pvalue[j,12] = sum(m22)/200
  org_result_pvalue[j,13] = sum(m23)/200
  
  org_result_pvalue[j,14] = sum(m31)/200
  org_result_pvalue[j,15] = sum(m32)/200
  org_result_pvalue[j,16] = sum(m33)/200
  
  
  org_result_pvalue[j,17] = sd(m11)/sqrt(200)
  org_result_pvalue[j,18] = sd(m12)/sqrt(200)
  org_result_pvalue[j,19] = sd(m13)/sqrt(200)
  
  org_result_pvalue[j,20] = sd(m21)/sqrt(200)
  org_result_pvalue[j,21] = sd(m22)/sqrt(200)
  org_result_pvalue[j,22] = sd(m23)/sqrt(200)
  
  org_result_pvalue[j,23] = sd(m31)/sqrt(200)
  org_result_pvalue[j,24] = sd(m32)/sqrt(200)
  org_result_pvalue[j,25] = sd(m33)/sqrt(200)
}

org_result_pvalue

