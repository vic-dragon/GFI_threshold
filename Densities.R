#######################
# Parameter Densities # 
#######################

#' new.threshold ~ Truncated N(old.threshold, small variance) 
#' old.threshold ~ Truncated N(new.threshold, small variance)
#' Alphas ~ N((XtX)^-1Xty, (XtX)^-1)
#' sigma^s ~ Inv-Chisquare(Scale=s2, df=n-3)   

## Distribution of Threshold

rnorm.trunc <- function(n, mean, sd, min, max){
  
  # Generate Random number from the truncated normal distribution
  
  r.num <- runif(n, min=pnorm(min, mean, sd), max=pnorm(max, mean, sd))
  res <- qnorm(r.num, mean, sd)
  
  return(res)
}

dnorm.trunc <- function(x, mean, sd, min, max){
  
  # Return log-probability value from the truncated normal distribution
  
  res<-dnorm(x, mean, sd) / (pnorm(max, mean, sd) - pnorm(min, mean, sd))
  
  return(log(res))
}

## Distribution of sigma

rchisq.inv.scaled <- function(n, scale, df){
  
  # Generate Random number from the scaled inverse chisquare dist.
  
  res <- scale*df/rchisq(n, df)
  
  return(res)
}

dchisq.inv.scaled <- function(x, scale, df){
  
  # Return log-probability value from the scaled inverse chisquare dist.
  
  res <- 0.5*df*log(0.5*df) + 0.5*df*log(scale) - lgamma(0.5*df) - (0.5*df+1)*log(x) - 0.5*df*scale/x
  
  return(res)
}