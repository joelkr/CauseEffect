
# Squared Error function used to calculate covariance
k <- function(Xi,Xj, p) p[1]^2 * exp(-0.5 * abs((Xi - Xj)/ p[2])^2)

# Linear covariance
kl <- function(Xi, Xj, p) p[1]^2 + p[2] * (Xi * Xj)

# Could try using this:
# Parameters:
# X1, X2 = vectors
# l = the scale length parameter
# Returns:
# a covariance matrix
cov <- function(X1,X2, f=k, p) {
#cov <- function(X1,X2, v, l) {
#cov <- function(X1,X2, params) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      # do k function here. remove comments. he's profiled for speed.
      #Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
      Sigma[i,j] <- f(X1[i], X2[j], p)
      #Sigma[i,j] <- k(X1[i], X2[j], params)
    }
  }
  return(Sigma)
}

gpK <- function(params, varN, x) {
  #print("gpK x:")
  #print(x)
  #print("gpK params:")
  #print(params)
  sigF <- params['sigF']
  l2   <- params['l2']
  sigL <- params['sigL']
  l1   <- params['l1']
  a    <- params['a']
  b    <- params['b']
  
  # Calculate K^-1 (K inverse)
  #K <- (cov(x, x, sigF=params['sigF'], l=params['l'])) + varN * diag(1, length(x))
  #K <- (cov(x, x, params)) + varN * diag(1, length(x))
  K <- (cov(x, x, f=k, c(sigL, l1))) + (cov(x, x, f=k, c(sigF, l2))) + varN * diag(1, length(x))#+ cov(x, x, f=kl, c(b, a)) + varN * diag(1, length(x))
  #Kinv <- solve(K)

  #return(list(K=K, Kinv=Kinv))
  return(K)
}

# all the data: x, y, varN,llC may need to be passed in as dataframe. 
# Kluge: tack all the scalars down the last column and keep track of which is 
# which. obs$p[1] == varN, obs$p[2] == llC
maxLogLik <- function(params, obs) {
  #print("maxLogLik params:")
  #print(params)
  #print("******maxLogLik obs*******:")
  #print(obs)

  logLikY <- function(params, obs)  {
    #with(obsa, { # forgot brace
    #Klist <- gpK(params=params, varN=obs$p[1], x=obs$x)
    K <- gpK(params=params, varN=obs$p[1], x=obs$x)
    #K <- cov(params, obs$x)
    #z <- obsa$y * 2
    # This should be: y %*% solve(K,y)
    #llY <- -(0.5 * t(obs$y) %*% Klist$Kinv %*% obs$y) - (0.5 * log(det(Klist$K))) - obs$p[2]
    llY <- -(0.5 * obs$y %*% solve(K, obs$y)) - (0.5 * log(abs(det(K)))) - obs$p[2]
    # })
    return(-llY)
  }

  return(optim(par=params, logLikY, obs=obs)$par)
  #return(0)
}

gpPredictEf <- function(params, obs, x_predict)  {
  #Klist <- gpK(params=params, varN=obs$p[1], x=obs$x)
  K <- gpK(params=params, varN=obs$p[1], x=obs$x)
  # ybar = Kstar * Kinv * y # best estimate of y is mean of distribution
  #Kstar <- cov(x_predict, obs$x, params['sigF'], params['l'])
  Kstar <- cov(x_predict, obs$x, f=k, c(params['sigL'], params['l1'])) + cov(x_predict, f=k, obs$x, c(params['sigF'], params['l2'])) #+ cov(x_predict, obs$x, f=kl, c(params['b'], params['a']))
  #Ef <- cov(x_predict, obs$x, obs$p[2], params['l']) %*% Klist$Kinv %*% obs$y
  # t(y) %*% solve(K) %*% y should be done: y %*% solve(K,y)
  #Ef <- Kstar %*% Klist$Kinv %*% obs$y
  Ef <- Kstar %*% solve(K, obs$y)

  return(Ef)
}

setupObs <- function(x, y, varN=0) {
  n <- length(x)
  llC <- (n/2) * log(2*pi)
  obs <- data.frame(x = x, y = y, p = rep(0, n))
  # Experiment did not provide error estimate.
  if(varN == 0)  {
  #varN <- 0.3 ^ 2
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

#setupParams <- function(x, y) {
setupParams <- function(x, y, a = 0, b=0) {
  # k needs to have and estimate for l which seems to be standard deviation
  # of the sampling grid.
  l1 <- sd(diff(sort(x))) # Local width of gaussian
  if(l1 == 0) l1 <- diff(sort(x))[1]/1000
  sigL <- sd(diff(sort(y))) # Local height of gaussian
  if(sigL == 0) sigL <- diff(sort(y))[1]
  sigF <- sd(y) # Global height of gaussian
  l2 <- sd(x) # Global width of gaussian
#params$l <- l
#params <- c(sigF = sigF, l = l, varN = varN)
  #params <- c(sigF = sigF, l1 = l1, sigL = sigL, l2 = l2)
  params <- c(sigF = sigF, l1 = l1, sigL = sigL, l2 = l2, a=a, b=b)
  return(params)
}

#x <- c(-1.5, -1.0, -0.75, -0.40, -0.25, 0.0)
#y <- c(-1.4, -0.9, -0.5,   0.2,   0.5,  0.9)
#varN <- 0.3
#x <- pressure$temperature
#y <- pressure$pressure
#varN <- 0
x <- trees$Height
y <- trees$Girth
varN <- 0
# x_predict is the sampling grid
#x_predict <- seq(min(x), max(x), len=(10*length(x))) # Probably not right
x_predict <- seq(min(x), max(x), len=(length(x))) # Probably not right
x_predict <- x
obs <- setupObs(x, y, varN)
# Get lms slope and intercept
lm1 <- lm(obs$y ~ obs$x)
a <- lm1$coeff[2]
b <- lm1$coeff[1]
names(a) <- ""
names(b) <- ""
params <- setupParams(x, y, a=a, b = b)
# Run a test and plot with the tutorial data.
print("Estimating parameters.")
p2 <- maxLogLik(params, obs)
print("Calculating expected value function.")
Ef <- gpPredictEf(p2, obs, x_predict)

plot(obs$x, obs$y, col='red')
points(x_predict, Ef, col='green')
lines(x_predict, Ef, col='green')
lines(obs$x, lm1$fitted, col='blue')

