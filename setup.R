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
params <- 0

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
