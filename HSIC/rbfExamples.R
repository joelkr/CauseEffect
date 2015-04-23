rbf.gauss <- function(gamma=1.0) {

  function(x) {
    exp(-gamma * norm(as.matrix(x),"F") )
  }
}

D <- matrix(c(-3,1,4), ncol=1) # 3 datapoints
N <- length(D[,1])

xlim  <- c(-5,7)
plot(NULL,xlim=xlim,ylim=c(0,1.25),type="n")
points(D,rep(0,length(D)),col=1:N,pch=19)
x.coord = seq(-7,7,length=250)
gamma <- 1.5
for (i in 1:N) {
  points(x.coord, lapply(x.coord - D[i,], rbf.gauss(gamma)), type="l", col=i)
}

# returns a rbf model given the:
# * observations x1, xN of dataset D
# * output value for each observation
# * number of centers
# * gamma value

rbf <- function(X, Y, K=10, gamma=1.0) {
  N     <- dim(X)[1] # number of observations
  ncols <- dim(X)[2] # number of variables

  repeat {
    km <- kmeans(X, K)  # let's cluster K centers out of the dataset
    if (min(km$size)>0) # only accept if there are no empty clusters
      break
  }

  mus <- km$centers # the clusters points

  Phi <- matrix(rep(NA,(K+1)*N), ncol=K+1)
  for (lin in 1:N) {
    Phi[lin,1] <- 1    # bias column
    for (col in 1:K) {
      Phi[lin,col+1] <- exp( -gamma * norm(as.matrix(X[lin,]-mus[col,]),"F")^2 )
    }
  }

  w <- pseudoinverse(t(Phi) %*% Phi) %*% t(Phi) %*% Y  # find RBF weights

  list(weights=w, centers=mus, gamma=gamma)  # return the rbf model
}

# An implementation for the prediction function

rbf.predict <- function(model, X, classification=FALSE) {

  gamma   <- model$gamma
  centers <- model$centers
  w       <- model$weights
  N       <- dim(X)[1]    # number of observations

  pred <- rep(w[1],N)  # we need to init to a value, so let's start with the bias

  for (j in 1:N) {  
    # find prediction for point xj
    for (k in 1:length(centers[,1])) {
      # the weight for center[k] is given by w[k+1] (because w[1] is the bias)
      pred[j] <- pred[j] + w[k+1] * exp( -gamma * norm(as.matrix(X[j,]-centers[k,]),"F")^2 )
    }
  }

  if (classification) {
    pred <- unlist(lapply(pred, sign))
  }
  pred
}


