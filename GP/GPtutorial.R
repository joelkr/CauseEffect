source("GPToy.R")

#x <- c(-1.5, -1.0, -0.75, -0.40, -0.25, 0.0)
#y <- c(-1.4, -0.9, -0.5,   0.2,   0.5,  0.9)
#varN <- 0.3
#x <- pressure$temperature
#y <- pressure$pressure
#varN <- 0
#x <- trees$Height
#y <- trees$Girth
x <- Indometh$time
y <- Indometh$conc
varN <- 0
# x_predict is the sampling grid
#x_predict <- seq(min(x), max(x), len=(10*length(x))) # Probably not right
obs <- setupObs(x, y, varN)
#x_predict <- seq(min(x), max(x), len=(length(x))) # Probably not right
# obs is sorted on x in setupObs()
# Get lms slope and intercept
lm1 <- lm(obs$y ~ obs$x)
a <- lm1$coeff[2]
b <- lm1$coeff[1]
names(a) <- ""
names(b) <- ""
#params <- setupParams(x, y, a=a, b = b)
params <- setupParams(x, y)
# Run a test and plot with the tutorial data.
print("Estimating parameters.")
p2 <- maxLogLik(params, obs)
#x_predict <- obs$x
#x_predict <- undiscretizeX(obs$x, params['lL'])
x_predict <- seq(min(obs$x), max(obs$x), len=length(obs$x)*10)
print("Calculating expected value function.")
Ef <- gpPredictEf(p2, obs, x_predict)

plot(obs$x, obs$y, col='red')
points(x_predict, Ef, col='green')
lines(x_predict, Ef, col='green')
lines(obs$x, lm1$fitted, col='blue')

