Fid.density <- function(x, y, threshold, alphas, sigma2, order, lambda){
  
  # Fiducial Density for the Spline model
  # num.case 1 : Constant Spline
  # num.case 2 : Linear Spline
  # num.case 3 : Quadratic Spline
  
  s=sqrt(sigma2)
  
  elambda<-exp(-lambda*(x-threshold))
  
  Col2<-(x-threshold)^order/(1+elambda)
  
  X <- cbind(rep(1,length(x)),Col2) # Design Matrix 
  
  l <- sum(dnorm((y - X %*% alphas),sd=s,log=T)) # Likelihood
  
  J <- Cal.J(x, y, alphas, threshold, s, order, lambda) # Determinant of Jacobian

  if(sum(is.na(J))>0){
    J <- Cal.J.indicator(x, y, alphas, threshold, s, order)  
  }
  
  res <- l + log(J)
  
  return(res)
}



Cal.J<-function(x, y, alpha.hat, threshold, s, order, lambda){
  
  #' Calculate Determinant of Jacobian
  #' Since it takes much times to consider all the tuples, 
  #' we will consider only 100 tuples and calculate the determinant for it

  elambda<-exp(-lambda*(x-threshold))
  
  Col2<-(x-threshold)^order/(1+elambda)
  
  Col3<-order*(x-threshold)^(order-1)/(1+elambda)*(order-1>=0) + lambda*elambda/(1+elambda)^2*(x-threshold)^order
  
  J<-cbind(rep(1,length(x)), Col2, Col3, y)
  
  tuples <- sapply(1:100, function(x) sample(dim(J)[1],dim(J)[2]))
  det.jaco <- sapply(1:100, function(x) det(J[tuples[,x],]))
  
  res <- sum(abs(det.jaco))*as.vector(alpha.hat[2])/s
  
  return(abs(res))
}

Cal.J.indicator<-function(x, y, alpha.hat, threshold, s, order){
  
  #' Calculate Determinant of Jacobian
  #' Since it takes much times to consider all the tuples, 
  #' we will consider only 100 tuples and calculate the determinant for it
  
  if(order==0){
    J <- cbind(rep(1,length(x)),(x>=threshold),y)
    # Jacobian Matrix for Constant Spline
  }else{
    J <- cbind(rep(1,length(x)),(x-threshold)^(order)*(x>=threshold),(x-threshold)^(order-1)*(x>=threshold),y)
  }
  
  tuples <- sapply(1:100, function(x) sample(dim(J)[1],dim(J)[2]))
  det.jaco <- sapply(1:100, function(x) det(J[tuples[,x],]))
  
  res <- sum(abs(det.jaco))*as.vector(alpha.hat[2])/s
  
  return(abs(res))
}



#############################
# Parameters for Prior dist #
#############################
prior.info <- function(x, y, threshold, df, order, lambda){
  
  # Calculate the parameters for the prior distribution

  elambda <- exp(-lambda*(x-threshold))
  Col2 <- (x-threshold)^order/(1+elambda)
  X <- cbind(rep(1,length(x)),Col2) # Design Matrix 
  
  alpha.var <- chol2inv(chol(crossprod(X)))
  svd.alpha.var<-svd(alpha.var)
  alpha.sd<-svd.alpha.var$u %*% diag(sqrt(svd.alpha.var$d)) %*% t(svd.alpha.var$v)
  alpha.hat <- alpha.var %*% t(X) %*% y
  
  s2 <- crossprod(y - X %*% alpha.hat)/df
  
  return(list(alpha.hat=alpha.hat, alpha.sd=alpha.sd, s2=as.vector(s2)))
}

###########################
# Metropolis-Hasting alg. # 
###########################

Metro.alg.spline<-function(x, y, lambda, old.threshold=0.5, sd.threshold = 0.1, burn.in=500, s.size=1000){
  
  res.store<-matrix(ncol=5,nrow=s.size)
  
  x.range=range(x)
  
  df <- length(y) - 3 # the number of parameter is 3 (intercept, 1 coefficient, 1 threshold)
  
  prob.model=f.dense.model.selection(x, y, threshold=old.threshold, lambda=lambda)
  selected.model=sample(1:4,1,prob=unlist(prob.model))
  
  # Start point
  old.prior.info <- prior.info(x=x, y=y, threshold=old.threshold, df=df, order=selected.model-1, lambda=lambda)
  
  old.var <- rchisq.inv.scaled(1,old.prior.info$s2,df)
  old.z <- rnorm(2)
  old.alpha <- old.prior.info$alpha.hat + sqrt(old.prior.info$s2) * old.prior.info$alpha.sd %*% old.z
  
  old.fid  <- Fid.density(x, y, old.threshold, old.alpha, old.var, order=selected.model-1, lambda=lambda)
  
  temp.alpha0<-c()
  temp.alpha1<-c()
  temp.var<-c()
  temp.threshold<-c()
  temp.model<-c()
  
  
  for( i in 1:(burn.in+s.size) ){
    
    new.threshold <- rnorm.trunc(1, mean=old.threshold, sd=sd.threshold, min=x.range[1], max=x.range[2])

    new.prob.model=f.dense.model.selection(x,y,threshold=new.threshold, lambda=lambda)
    new.selected.model=sample(1:4,1,prob=unlist(new.prob.model))
    
    
    new.prior.info <- prior.info(x=x, y=y, threshold=new.threshold, df=df, order=new.selected.model-1, lambda=lambda)
    
    new.var<-rchisq.inv.scaled(1,new.prior.info$s2,df)
    new.z<- rnorm(2)
    new.alpha<- new.prior.info$alpha.hat + sqrt(new.prior.info$s2) * new.prior.info$alpha.sd %*% new.z
    
    new.fid <- Fid.density(x, y, new.threshold, new.alpha, new.var, order=new.selected.model-1, lambda=lambda)
    
    new.prop.dense <-  dnorm.trunc(new.threshold, mean=old.threshold, sd=sd.threshold, min=x.range[1], max=x.range[2])+
      sum(dnorm(new.z, log=TRUE)) + dchisq.inv.scaled(new.var, new.prior.info$s2, df)
    
    old.prop.dense <-  dnorm.trunc(old.threshold, mean=new.threshold, sd=sd.threshold, min=x.range[1], max=x.range[2])+
      sum(dnorm(old.z, log=TRUE)) + dchisq.inv.scaled(old.var, old.prior.info$s2, df)
    
    r <- exp( new.fid - old.fid + old.prop.dense - new.prop.dense)
    
    if( runif(1) < r ) {
      
      temp.threshold<-c(temp.threshold, new.threshold)
      
      temp.alpha0<-c(temp.alpha0, as.vector(new.alpha)[1])
      temp.alpha1<-c(temp.alpha1, as.vector(new.alpha)[2])
      
      temp.var<-c(temp.var, new.var)
      temp.model<-c(temp.model, new.selected.model)
      
      old.fid <- new.fid
      
      old.threshold <- new.threshold
      
      old.alpha <- new.alpha
      
      old.var <- new.var
      
      old.z<-new.z
      
      old.prior.info<-new.prior.info
      
      selected.model<-new.selected.model
      
    }else{
      temp.threshold<-c(temp.threshold, old.threshold)
      
      temp.alpha0<-c(temp.alpha0, as.vector(old.alpha)[1])
      temp.alpha1<-c(temp.alpha1, as.vector(old.alpha)[2])
      
      temp.var<-c(temp.var, old.var)
      temp.model<-c(temp.model, selected.model)
    }
  }
  res.store[,1]<-temp.alpha0[(burn.in+1):(burn.in+s.size)]
  res.store[,2]<-temp.alpha1[(burn.in+1):(burn.in+s.size)]
  res.store[,3]<-temp.var[(burn.in+1):(burn.in+s.size)]
  res.store[,4]<-temp.threshold[(burn.in+1):(burn.in+s.size)]
  res.store[,5]<-temp.model[(burn.in+1):(burn.in+s.size)]
  
  return(res.store)
}
