# Squared Error function used to calculate covariance
k <- function(Xi,Xj, p) p[1]^2 * exp(-0.5 * abs((Xi - Xj)/ p[2])^2)

# Partials for parameter 1 and 2
p1k <- function(Xi, Xj, p) 2 * p[1] * exp(-0.5 * abs((Xi - Xj)/p[2])^2)
p2k <- function(Xi, Xj, p) p[1] ^ 2 * exp((Xi - Xj)^2/p[2]^2)*((Xi - Xj)^2/p[2]^3)

# Linear covariance
kl <- function(Xi, Xj, p) p[1]^2 + p[2] * (Xi * Xj)

cov <- function(X1,X2, f=k, p) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- f(X1[i], X2[j], p)
    }
  }
  return(Sigma)
}

# K might be calculated:
# as.matrix(dist(x, method="euclidean", diag=TRUE))
# This would be new cov() function
# a^2 - 2*a*b + b^2 is less stable than (a-b)^2 because of numerical precision.
# This means one should subtract the mean first. This shouldn't change anything
# because squared error is independent of the mean.
# May need to create a different covariance function for each k function.
cov2 <- function(X1, X2, f=k, p) {
  ######
  ###### Need to fix this so that it will do the right thing when X1 and X2
  # are different. dist() treats every row in a matrix as a vector and returns
  # the distance between every pair of vectors.
  x <- X1 - mean(X1)
  K <- p[1]^2 * exp(-0.5 * as.matrix(dist(x, method="euclidean", diag=TRUE))^2/p[2]^2)
  return(K)
}

gpK <- function(params, varN, x) {
  sigF <- params['sigF']
  lF   <- params['lF']
  sigL <- params['sigL']
  lL   <- params['lL']
#  varN <- params['varN']
  # Remove a and b if linear covariance is not going to be used.
  #a    <- params['a']
  #b    <- params['b']
  
  #K <- (cov(x, x, f=k, c(sigL, lL))) + (cov(x, x, f=k, c(sigF, lF))) + varN * diag(1, length(x))#+ cov(x, x, f=kl, c(b, a))
  K <- (cov2(x, x, f=k, c(sigL, lL))) + (cov2(x, x, f=k, c(sigF, lF))) + varN * diag(1, length(x))#+ cov(x, x, f=kl, c(b, a))

  return(K)
}

gpGrad <- function(p_i, i, varN, x) {
  if(i == 1) return(cov(x,x, f=p1k, p_i) + varN * diag(1, length(x)))
  else if(i == 2) return(cov(x,x, f=p2k, p_i) + varN * diag(1, length(x)))
  else return(NA)
}
# all the data: x, y, varN,llC may need to be passed in as dataframe. 
# Kluge: tack all the scalars down the last column and keep track of which is 
# which. obs$p[1] == varN, obs$p[2] == llC
# These may be right. Then need to generate each derivative K matrix.
# partial k/p1 = 2 * p1 * exp(-0.5(Xi - Xj)^2/p2^)
# partial k/p2 = p1^2 * exp((Xi - Xj)^2/p2)
maxLogLik <- function(params, obs) {

  logLikY <- function(params, obs)  {
    K <- gpK(params=params, varN=obs$p[1], x=obs$x)
    llY <- -(0.5 * obs$y %*% solve(K, obs$y)) 
             - (0.5 * log(abs(det(K)))) - obs$p[2]
    return(-llY)
  }
  grLogLikY <- function(params, obs) {
    grad <- rep(0, length(params))

    return(grad)
  }
  return(optim(par=params, logLikY, obs=obs)$par)
}

gpPredictEf <- function(params, obs, x_predict)  {
  K <- gpK(params=params, varN=obs$p[1], x=obs$x)
  #Kstar <- cov(x_predict, obs$x, f=k, c(params['sigL'], params['lL'])) + cov(x_predict, f=k, obs$x, c(params['sigF'], params['lF'])) #+ cov(x_predict, obs$x, f=kl, c(params['b'], params['a']))
  Kstar <- cov2(x_predict, obs$x, f=k, c(params['sigL'], params['lL'])) + cov2(x_predict, f=k, obs$x, c(params['sigF'], params['lF'])) #+ cov(x_predict, obs$x, f=kl, c(params['b'], params['a']))
  # Formula: ExpectedVal <- Kstar %*% K^-1 %*% y
  Ef <- Kstar %*% solve(K, obs$y)

  return(Ef)
}

undiscretizeX <- function(x, lL) {
  x <- sort(x)
  s <- lL/100
  z <- rep(0, length(x))
  i <- 2
  while(i <= length(x)) {
    if(x[i] == x[i-1]) {
      z[i] <- z[i] + rnorm(1, 0, s)
    }
    i <- i+1
  }
  return(x + z)
}

setupObs <- function(x, y, varN=0) {
#setupObs <- function(x, y) {
  n <- length(x)
  llC <- (n/2) * log(2*pi)
  obs <- data.frame(x = x, y = y, p = rep(0, n))
  # Sort by presumed independent variable, may want a flag to control.
  obs <- obs[order(obs[,1]), ]
  # Experiment did not provide error estimate.
  if(varN == 0)  {
    varN <- (var(diff(sort(x))/100))
    # Experiments sample x or t at an even interval.
    if(varN == 0) {
      varN <- (diff(sort(x))[1]/100)^2
    }
  }
  # Must do this AFTER sorting and scaling
  obs$p[1] <- varN
  obs$p[2] <- llC
  return(obs)
}

#setupParams <- function(x, y, a = 0, b=0) {
setupParams <- function(x, y) {
#setupParams <- function(x, y, varN=0) {
  # k needs to have and estimate for l which seems to be standard deviation
  # of the sampling grid.
  lL <- sd(diff(sort(x))) # Local width of gaussian
  # Hopefully less than error in sampling grid_x. Must add some noise or
  # many K matrices seem to be impossible to solve.
  if(lL == 0) lL <- diff(sort(x))[1]/100
  sigL <- sd(diff(sort(y))) # Local height of gaussian
  if(sigL == 0) sigL <- diff(sort(y))[1]
  sigF <- sd(y) # Global height of gaussian
  lF <- sd(x) # Global width of gaussian
  # Experiment did not provide error estimate.
#  if(varN == 0)  {
#    varN <- (var(diff(sort(x))/100))
#    # Experiments sample x or t at an even interval so var is zero.
#    if(varN == 0) {
#      # Arbitrarily decide they can accurately measure 1/100th of the 
#      # interval.
#      varN <- (diff(sort(x))[1]/100)^2
#    }
#  }
  #params <- c(sigF = sigF, lL = lL, sigL = sigL, lF = lF, a=a, b=b)
  params <- c(sigL = sigL, lL = lL, sigF = sigF, lF = lF)
#  params <- c(sigL = sigL, lL = lL, sigF = sigF, lF = lF, varN=varN)
  return(params)
}

