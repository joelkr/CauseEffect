# This file should try using a simple linear model for regression
source("hsic.R")


pd <- read.table("EData/pair0001.txt")
obs <- data.frame(x = pd$V1, y = pd$V2, xs = scale(pd$V1), ys = scale(pd$V2))

# First we want to try doing things without sorting or scaling
X <- matrix(obs$x, length(obs$x), 1)
Y <- matrix(obs$y, length(obs$y), 1)
Xscaled <- matrix(obs$xs, length(obs$x), 1)
Yscaled <- matrix(obs$ys, length(obs$y), 1)

# Linear model of X -> Y
lmXtoY <- lm(obs$y ~ obs$x)
# Linear model of Y -> X
lmYtoX <- lm(obs$x ~ obs$y)

# Linear model of X -> Y scaled
scaledLmXtoY <- lm(obs$ys ~ obs$xs)
# Linear model of Y -> X scaled
scaledLmYtoX <- lm(obs$xs ~ obs$ys)

# Residuals
eY <- matrix(lmXtoY$resid, length(lmXtoY$resid), 1)
eX <- matrix(lmYtoX$resid, length(lmYtoX$resid), 1)

eYscaled <- matrix(scaledLmXtoY$resid, length(scaledLmXtoY$resid), 1)
eXscaled <- matrix(scaledLmYtoX$resid, length(scaledLmYtoX$resid), 1)

X_Y <- hsic(X, eY)
Y_X <- hsic(Y, eX)

X_Yscaled <- hsic(Xscaled, eYscaled)
Y_Xscaled <- hsic(Yscaled, eXscaled)

