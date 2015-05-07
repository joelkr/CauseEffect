# Translated from MATLAB code written by:
# Joris Mooij  <j.m.mooij@uva.nl>
# https://staff.fnwi.uva.nl/j.m.mooij/
require(pracma)

hsic <- function(X, Y, sX=0, sY=0) {
  N <- dim(X)[1]
  if(dim(X)[2] > 1 || dim(Y)[2] > 1) {
    print("Only column vectors for now")
    return(0)
  }
  if(dim(Y)[1] != N) {
    print("X and Y must have same number of rows")
    return(0)
  }
  if(sX == 0) {
    sX <- guess_sigma(X)
  }
  if(sY == 0) {
    sY <- guess_sigma(Y)
  }
  KX <- rbfkernel(X, sX)
  KY <- rbfkernel(Y, sY)

  H <- eye(N) - (1/N) * ones(N)
  KXbar <- H %*% KX %*% H
  KYbar <- H %*% KY %*% H

  HSIC <- 1/(N^2) * sum(diag(KXbar %*% KY))

  # calculate sums of kernel matrices
  KX_sum <- sum(KX)
  KY_sum <- sum(KY)

  # calculate statistics for gamma approximation
  x_mu <- 1.0 / (N * (N-1)) * (KX_sum - N)
  y_mu <- 1.0 / (N * (N-1)) * (KY_sum - N)
  mean_H0 <- (1.0 + x_mu * y_mu - x_mu - y_mu) / N
  var_H0 <- (2.0 * (N-4) * (N-5)) / (N * (N-1.0) * (N-2) * (N-3) * ((N-1)^4)) * sum(diag(KXbar %*% KX)) * sum(diag(KYbar %*% KY))

  #cat("x_mu: ", x_mu, "y_mu: ", y_mu, "mean_H0: ", mean_H0, "var_H0: ", var_H0, "\n")

  # calculate p-value under gamma approximation
  a <- mean_H0 * mean_H0 /var_H0
  b <- N * var_H0 / mean_H0
  pHSIC <- pgamma(N * HSIC / b, a, lower=FALSE)
  cat("N: ", N, "a: ", a, "b: ", b, "pHSIC: ", pHSIC, "\n")

  rv <- list()
  rv[["pHSIC"]] <- pHSIC
  rv[["HSIC"]] <- HSIC
  return(rv)

}

rbfkernel <- function(X, sigma) {
  # Make this work for vectors only for the moment.
  N <- dim(X)[1]
  d <- dim(X)[2]

  if(d > 1) {
    print("Only column vectors for now")
    return(0)
  }
  else {
    K <- (repmat(X,1,N) - repmat(t(X),N,1))^2
  }
  K <- exp(-K / (2.0 * sigma^2))
  return(K)
}

guess_sigma <- function(X, method=0) {
  if (method == 0) {
    Xnorm <- get_norm(X)
    # ?
    #Xnorm <- Xnorm - tril(Xnorm)
    # Try this:
    Xnorm[lower.tri(Xnorm, diag = TRUE)] <- 0
    #Xnorm <- matrix(Xnorm, dim(X)[1]^2, 1)
    #sigma <- sqrt(0.5 * median(Xnorm[Xnorm > 0]))
    # And this ( skip the reshaping)
    sigma <- sqrt(0.5 * median(Xnorm[Xnorm > 0]))
  }
  else if(method == 1){
    Xnorm <- get_norm(X)
    sigma <- sqrt(0.5*median(Xnorm))
  }
  else if(method == 2) {
    #sigma = exp(fminsearch(@(logh) kernel_LOO(exp(logh),X),0.0));
    print("Not implemented yet.")
    sigma <- NA
  }
  return(sigma)
}

get_norm <- function(A) {
  # A must be promoted to matrix
  lenA <- dim(A)[1]
  result <- zeros(lenA, lenA)
    for(i1 in 1:(lenA-1)) {
      for(i2 in (i1+1):lenA) {
        result[i1, i2] <- sum((A[i1,] - A[i2,])^2)
        #cat("result[", i1, ",", i2, "]: ", result[i1, i2], "\n")
        result[i2, i1] <- result[i1, i2]
      }
    }
  return(result)
}
