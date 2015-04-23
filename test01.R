print("X -> Y")
obs <- setupObs(edata[[1]][,1], edata[[1]][,2])
params <- setupParams(obs$x, obs$y)
print("Calculating optimal params")
system.time(p2 <-  maxLogLik(params, obs))
x_predict <- obs$x
print("Calculating expected Y")
EfX <- gpPredictEf(params, obs, x_predict)
# Residuals with X as cause
eX <- obs$y - EfX
print("Y -> X")
obs <- setupObs(edata[[1]][,2], edata[[1]][,1])
params <- setupParams(obs$x, obs$y)
print("Calculating optimal params")
system.time(p2 <-  maxLogLik(params, obs))
x_predict <- obs$x
print("Calculating expected Y")
EfY <- gpPredictEf(params, obs, x_predict)
# Residuals with X as cause
eY <- obs$y - EfY

