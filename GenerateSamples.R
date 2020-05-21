#Functions
## mu functions (non-decreasing)
f1<-function(x){
  (x>=0.5)*0.5
}
f2<-function(x){
  2*(x>=0.5)*(x-0.5)^2
}
f3<-function(x){
  lambda=0.5*log(2)
  (x>=0.5)*exp(-lambda/(x-0.45))
}

####################
# Generate Samples # 
####################

gen.sample<-function(n, m, fnct, sigma){
  #' n : the number of dose level
  #' m : the number of replication
  
  i <- 1:n
  x <- i/(1+n) #Dose level
  
  res.sample <- matrix(ncol=n,nrow=m)
  
  for(j in 1:n){
    exp.eval <- paste0("f",fnct,"(","x[",j,"]",")")
    res.sample[,j] <- eval(parse(text=exp.eval))+rnorm(m, mean=0, sd=sigma)
  }
  
  res <- cbind(rep(x,each=m),as.vector(res.sample))
  
  return(res)
}