# Squared Error function used to calculate covariance
k <- function(Xi,Xj, p) p[1]^2 * exp(-0.5 * abs((Xi - Xj)/ p[2])^2)

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

gpK <- function(params, varN, x) {
  sigF <- params['sigF']
  l2   <- params['l2']
  sigL <- params['sigL']
  l1   <- params['l1']
  # Remove a and b if linear covariance is not going to be used.
  a    <- params['a']
  b    <- params['b']
  
  # Calculate K^-1 (K inverse)
  K <- (cov(x, x, f=k, c(sigL, l1))) + (cov(x, x, f=k, c(sigF, l2))) + varN * diag(1, length(x))#+ cov(x, x, f=kl, c(b, a))

  return(K)
}

# all the data: x, y, varN,llC may need to be passed in as dataframe. 
# Kluge: tack all the scalars down the last column and keep track of which is 
# which. obs$p[1] == varN, obs$p[2] == llC
maxLogLik <- function(params, obs) {

  logLikY <- function(params, obs)  {
    K <- gpK(params=params, varN=obs$p[1], x=obs$x)
    llY <- -(0.5 * obs$y %*% solve(K, obs$y)) - (0.5 * log(abs(det(K)))) - obs$p[2]
    return(-llY)
  }
  return(optim(par=params, logLikY, obs=obs)$par)
}

gpPredictEf <- function(params, obs, x_predict)  {
  K <- gpK(params=params, varN=obs$p[1], x=obs$x)
  Kstar <- cov(x_predict, obs$x, f=k, c(params['sigL'], params['l1'])) + cov(x_predict, f=k, obs$x, c(params['sigF'], params['l2'])) #+ cov(x_predict, obs$x, f=kl, c(params['b'], params['a']))
  Ef <- Kstar %*% solve(K, obs$y)

  return(Ef)
}

setupObs <- function(x, y, varN=0) {
  n <- length(x)
  llC <- (n/2) * log(2*pi)
  obs <- data.frame(x = x, y = y, p = rep(0, n))
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

setupParams <- function(x, y, a = 0, b=0) {
  # k needs to have and estimate for l which seems to be standard deviation
  # of the sampling grid.
  l1 <- sd(diff(sort(x))) # Local width of gaussian
  # Hopefully less than error in sampling grid_x. Must add some noise or
  # many K matrices seem to be impossible to solve.
  if(l1 == 0) l1 <- diff(sort(x))[1]/100
  sigL <- sd(diff(sort(y))) # Local height of gaussian
  if(sigL == 0) sigL <- diff(sort(y))[1]
  sigF <- sd(y) # Global height of gaussian
  l2 <- sd(x) # Global width of gaussian
  params <- c(sigF = sigF, l1 = l1, sigL = sigL, l2 = l2, a=a, b=b)
  return(params)
}

