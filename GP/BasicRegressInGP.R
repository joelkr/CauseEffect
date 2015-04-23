require(MASS)
require(reshape2)
require(ggplot2)

set.seed(12345)

x_predict <- seq(-5,5,len=50)
l <- 1

SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
cov <- function(X, Y) outer(X, Y, SE, l)

COV <- cov(x_predict, x_predict)
values <- mvrnorm(3, rep(0, length=length(x_predict)), COV)

dat <- data.frame(x=x_predict, t(values))
dat <- melt(dat, id="x")
head(dat)

fig2a <- ggplot(dat,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +   theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2a

#Posterior distribution given the data

obs <- data.frame(x = c(-4, -3, -1,  0,  2),
                  y = c(-2,  0,  1,  2, -1))

# Without noise
cov_xx_inv <- solve(cov(obs$x, obs$x))
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict)

# Take 3 random samples from posterior distribution
values <- mvrnorm(3, Ef, Cf)

dat <- data.frame(x=x_predict, t(values))
dat <- melt(dat, id="x")

fig2b <- ggplot(dat,aes(x=x,y=value)) +
  geom_ribbon(data=NULL, 
              aes(x=x_predict, y=Ef, ymin=(Ef-2*sqrt(diag(Cf))), ymax=(Ef+2*sqrt(diag(Cf)))),
              fill="grey80") +
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=NULL,aes(x=x_predict,y=Ef), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2b

# With noise
sigma.n <- 0.25
cov_xx_inv <- solve(cov(obs$x, obs$x) + sigma.n^2 * diag(1, length(obs$x)))
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict)

values <- mvrnorm(3, Ef, Cf)

dat <- data.frame(x=x_predict, t(values))
dat <- melt(dat, id="x")

fig2c <- ggplot(dat,aes(x=x,y=value)) +
  geom_ribbon(data=NULL, 
              aes(x=x_predict, y=Ef, ymin=(Ef-2*sqrt(diag(Cf))), ymax=(Ef+2*sqrt(diag(Cf)))),
              fill="grey80") + # Var
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=NULL,aes(x=x_predict,y=Ef), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2c

