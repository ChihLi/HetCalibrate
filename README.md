# HetCalibrate
This R package allows calibration parameter estimation for inexact computer models with heteroscedastic errors. More details can be found in [Sung, Barber, and Walker (2023+)](https://arxiv.org/abs/1910.11518).

You can install the package using `install_github` function as follows:
``` r
library(devtools)
install_github("ChihLi/HetCalibrate")
```
If you want to reproduce the results in [Sung, Barber, and Walker (2023+)](https://arxiv.org/abs/1910.11518), please go to the [reproducibilty webpage](https://github.com/ChihLi/HetCalibrate-Reproducibility).

An example is given below.
``` r
library(HetCalibrate)
set.seed(1)

##### setting #####
# computer model
f.sim <- function(x, cpara) {
 return(c(exp(x/10)*sin(x) - sqrt(cpara^2 - cpara + 1) * (sin(cpara*x)+cos(cpara*x))))
}
df.sim <- function(x, cpara) {
 return(c(-sqrt(cpara^2-cpara+1)*(x*cos(x*cpara)-x*sin(x*cpara))-((2*cpara-1)*(sin(x*cpara)+cos(x*cpara)))/(2*sqrt(cpara^2-cpara+1))))
}

# variance process - constant variance
var.f <- function(x) (0.01+0.2*(x-pi)^2)^2

# physical process
p.fun <- function(x) exp(x/10)*sin(x)

# true parameter
true.cpara <- optim(0, fn = function(g) {
 x.grid <- seq(0,2*pi,0.01)
 mean((p.fun(x.grid) - f.sim(x.grid, g))^2)
},
lower = -0.3, upper = 0.3, method = "L-BFGS-B")$par


# observed input
X0 <- seq(0,2*pi, length.out = 8)
# mean process
pmean <- p.fun(X0)
# variance process
var.y <- var.f(X0)
# number of replicates
n.rep <- rep(5, length(X0))

# setting for lower and upper bounds of parameters
cpara_min <- -0.3
cpara_max <- 0.3
cpara_init.vt <- c(-0.2, 0, 0.2)

# simulate X and Z
X <- matrix(rep(X0, n.rep), ncol = 1)
Z <- rep(0, sum(n.rep))
for(i in 1:length(X0)) {
 Z[(ifelse(i==1,0,sum(n.rep[1:(i-1)]))+1):sum(n.rep[1:i])] <- pmean[i] + rnorm(n.rep[i], 0, sd = sqrt(var.y[i]))
}

model <- mleHetCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
                        lower = 0.01*max(X0), upper = 2.5*max(X0),
                        init = list("cpara" = 0),
                        settings = list(checkHom = FALSE, linkThetas = "none"),
                        covtype = "Matern5_2", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)

print(cpara.Het.OGP <- model$cpara)
xgrid <- matrix(seq(min(X0), max(X0), length.out = 101), ncol = 1)
predictions.Het <- predict(x = xgrid, object =  model)

plot(X, Z, ylab = 'y', xlab = "x")
Z0 <- hetGP::find_reps(X, Z)$Z0
points(X0, Z0, pch = 20, cex = 1.2)
lines(xgrid, predictions.Het$mean, col = 'red', lwd = 2)
curve(p.fun, min(X0), max(X0), add = TRUE, col = 1, lty = 2, lwd = 1)
lines(xgrid, qnorm(0.025, predictions.Het$mean, sqrt(predictions.Het$sd2 + predictions.Het$nugs)),
     col = 3, lty = 3, lwd = 2)
lines(xgrid, qnorm(0.975, predictions.Het$mean, sqrt(predictions.Het$sd2 + predictions.Het$nugs)),
     col = 3, lty = 3, lwd = 2)
lines(xgrid, f.sim(xgrid, cpara.Het.OGP), col = 4, lty = 2, lwd = 2)
```

