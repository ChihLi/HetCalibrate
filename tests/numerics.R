library(HetCalibrate)
ns <- asNamespace('HetCalibrate')
helper <- function(name) get(name, ns)
set.seed(42)
x <- matrix(runif(12), 6, 2)
y <- matrix(runif(8), 4, 2)
h <- 1e-6
for (type in c('Gaussian', 'Matern5_2', 'Matern3_2')) {
  for (other in list(NULL, y)) {
    cv <- function(t) hetGP::cov_gen(x, other, theta = t, type = type)
    numerical <- (cv(0.7+h)-cv(0.7-h))/(2*h)
    analytic <- helper('partial_cov_gen')(x, theta=0.7, type=type,
                   arg='theta_k', X2=other) * cv(0.7)
    stopifnot(max(abs(numerical-analytic)) < 1e-7)
    numerical <- (cv(0.7*(1.2+h))-cv(0.7*(1.2-h)))/(2*h)
    analytic <- helper('partial_cov_gen')(x, theta=0.7, type=type,
                  arg='k_theta_g', k_theta_g=1.2, X2=other) * cv(0.7*1.2)
    stopifnot(max(abs(numerical-analytic)) < 1e-7)
  }
}
# Check matrix helpers against direct expressions, including one-row matrices.
for (n in c(1, 5)) {
 A <- matrix(rnorm(n*3), n, 3); B <- matrix(rnorm(3*n), 3, n)
 stopifnot(isTRUE(all.equal(helper('fast_diag')(A,B),diag(A%*%B))))
 D <- diag(n)
 stopifnot(isTRUE(all.equal(helper('add_diag')(D,rep(2,n)),3*diag(n))))
 stopifnot(identical(D,diag(n)))
}
stopifnot(isTRUE(all.equal(helper('fast_tUY2')(c(2,3,1),1:6),c(3,12,6))))
# Exercise every supported kernel with both fitting models.
f <- function(x,cpara) c(sin(x)+cpara)
df <- function(x,cpara) rep(1,nrow(x))
X <- matrix(rep(seq(0,2,length.out=6),each=4),ncol=1)
Z <- f(X,0.2)+rnorm(nrow(X),sd=0.1+0.05*X[,1])
for (type in c('Gaussian','Matern5_2','Matern3_2')) {
 for (fit in list(mleHomCalibrate,mleHetCalibrate)) {
  model <- fit(X,Z,cpara_min=-1,cpara_max=1,init=list(cpara=0),
     lower=0.1,upper=3,covtype=type,f.sim=f,df.sim=df,
     settings=list(checkHom=FALSE,linkThetas='none'))
  pred <- predict(model, matrix(c(0.2,0.9,1.7),ncol=1))
  stopifnot(all(is.finite(c(model$cpara,pred$mean,pred$sd2,pred$nugs))),
            all(pred$sd2 >= -1e-8), all(pred$nugs >= 0))
 }
}
