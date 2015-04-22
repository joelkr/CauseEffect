source("GP/GPToy.R")
source("HSIC/hsic.R")

files <- (Sys.glob("./EData/*[0-9].txt"))
edata <- lapply(files, function(x) read.table(x, header = FALSE))

# Variables we only need to calculate once
# Local variability
#lL <- 0
# Total variability
#fL <- 0
# Sigma_n = measurement noise
#sn <- 0
tEnv <- new.env()
# Grid values to predict x on
assign("x_predict", 0, envir = tEnv)
params <- list()

######
# Probably discretize x as well. Do this when preparing the dataframe
# originally. Scaling seems to make an uninvertible matrix. Probably want
# to have a separate setup file with the functions to set up each cause effect
# pair.
setupX <- function(x) {
  params <- list()
  # ###### Pull to separate function
  #        Make Global
  # Input equal length vectors x and y.
  # May need to add noise to x to make an invertible matrix.
  # Data needs to be sorted on x.
  # SE needs to have and estimate for l which seems to be standard deviation
  # of the sampling grid.
  # x must be sorted for this to work.
  x <- sort(x)
  #assign("lL", sd(diff(x)), envir = .GlobalEnv)
  params$lL <- sd(diff(x))
  # Second l for overall trend in the data
  params$fL <- sd(x)
  #assign("fL", (max(x) - min(x)), envir = .GlobalEnv)
  # To reliably invert x we need to add some small amount of possible error
  # to x. How much?
  # Another possibilty is to add a small random amount to every x that is 
  # equal - undiscretize.
  #assign("sn", (lL / 1000), envir = .GlobalEnv)
  params$sn <- params$lL/1000
  # x_predict is the sampling grid
  #x_predict <- seq(min(x), max(x), len=length(x)) # Probably not right
  assign("x_predict", seq(min(x), max(x), len=length(x)), envir = tEnv)
  # #######
  #
  return(params)
}

setupY <- function(y) {
  sigmaF <- sd(y)
  return(sigmaF)
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
