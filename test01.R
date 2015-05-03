source("setup.R")

print("X -> Y")
obsX <- setupObs(edata[[1]][,1], edata[[1]][,2])
params <- setupParams(obsX$x, obsX$y)
print("Calculating optimal params")
system.time(p2 <-  maxLogLik(params, obsX))
x_predict <- obsX$x
print("Calculating expected value with X as cause: EfX")
EfX <- gpPredictEf(p2, obsX, x_predict)
# Residuals with X as cause
print("Residuals with X as cause: eX")
eX <- obsX$y - EfX
print("Y -> X")
obsY <- setupObs(edata[[1]][,2], edata[[1]][,1])
params <- setupParams(obsY$x, obsY$y)
print("Calculating optimal params")
system.time(p2 <-  maxLogLik(params, obsY))
x_predict <- obsY$x
print("Calculating expected value with Y as cause: EfY")
EfY <- gpPredictEf(params, obsY, x_predict)
# Residuals with Y as cause
print("Residuals with Y as cause: eY")
eY <- obsY$y - EfY

