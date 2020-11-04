ddKg_ddphi_func <- function(object, arg = 1){
  out <- matrix(0, nrow(object$X0), nrow(object$X0))
  for(i in 1:nrow(object$X0)){
    for(j in i:nrow(object$X0)){
      x <- object$X0[i,arg,drop=TRUE]
      y <- object$X0[j,arg,drop=TRUE]
      theta <- object$theta_g[arg]
      if(object$covtype == "Matern5_2"){
        out[i,j] <- out[j,i] <-
          -(((y^2-2*x*y+x^2)*(15*theta^2+3*5^(3/2)*abs(y-x)*theta-25*y^2+50*x*y-25*x^2)*exp(-(sqrt(5)*abs(y-x))/theta))/(3*theta^6))
      }else if(object$covtype == "Matern3_2"){
        out[i,j] <- out[j,i] <-
          -((y-x)^2*(9*theta-3^(3/2)*abs(y-x))*exp(-(sqrt(3)*abs(y-x))/theta))/theta^5
      }else if(object$covtype == "Gaussian"){
        out[i,j] <- out[j,i] <-
          -((y-x)^2*(2*theta-y^2+2*x*y-x^2)*exp(-(x-y)^2/theta))/theta^4
      }
    }
  }
  return(out)
}

##'Compute Information matrix of a heteroscedastic model object (of class \code{hetCalibrate})
##' @param object an object of class \code{hetCalibrate}; e.g., as returned by \code{\link[HetCalibrate]{mleHetCalibrate}} with \code{orthogonal=TRUE}
##' @return a matrix which contains the information matrix of all the parameters. The order of the parameters is \code{cpara}, \code{theta},  \code{nu_hat}, \code{theta_g}, \code{g}, \code{nu_hat_var}, \code{Delta}
##' @examples
##' ##------------------------------------------------------------
##' ## Model calibration under heteroscedastic noises:
##' ##    model discrepancy is modeled by an orthogonal Gaussian process
##' ##------------------------------------------------------------
##'library(HetCalibrate)
##'set.seed(1)
##'
##'##### setting #####
##'# computer model
##'f.sim <- function(x, cpara) {
##'  return(c(exp(x/10)*sin(x) - sqrt(cpara^2 - cpara + 1) * (sin(cpara*x)+cos(cpara*x))))
##'}
##'df.sim <- function(x, cpara) {
##'  return(c(-sqrt(cpara^2-cpara+1)*(x*cos(x*cpara)-x*sin(x*cpara))-((2*cpara-1)*(sin(x*cpara)+cos(x*cpara)))/(2*sqrt(cpara^2-cpara+1))))
##'}
##'
##'# variance process - constant variance
##'var.f <- function(x) (0.01+0.2*(x-pi)^2)^2
##'
##'# physical process
##'p.fun <- function(x) exp(x/10)*sin(x)
##'
##'# true parameter
##'true.cpara <- optim(0, fn = function(g) {
##'  x.grid <- seq(0,2*pi,0.01)
##'  mean((p.fun(x.grid) - f.sim(x.grid, g))^2)
##'},
##'lower = -0.3, upper = 0.3, method = "L-BFGS-B")$par
##'
##'
##'# observed input
##'X0 <- seq(0,2*pi, length.out = 8)
##'# mean process
##'pmean <- p.fun(X0)
##'# variance process
##'var.y <- var.f(X0)
##'# number of replicates
##'n.rep <- rep(5, length(X0))
##'
##'# setting for lower and upper bounds of parameters
##'cpara_min <- -0.3
##'cpara_max <- 0.3
##'cpara_init.vt <- c(-0.2, 0, 0.2)
##'
##'# simulate X and Z
##'X <- matrix(rep(X0, n.rep), ncol = 1)
##'Z <- rep(0, sum(n.rep))
##'for(i in 1:length(X0)) {
##'  Z[(ifelse(i==1,0,sum(n.rep[1:(i-1)]))+1):sum(n.rep[1:i])] <- pmean[i] + rnorm(n.rep[i], 0, sd = sqrt(var.y[i]))
##'}
##'
##'model <- mleHetCalibrate(X = X, Z = Z, cpara_max = cpara_max, cpara_min = cpara_min,
##'                         lower = 0.01*max(X0), upper = 2.5*max(X0),
##'                         init = list("cpara" = 0),
##'                         settings = list(checkHom = FALSE, linkThetas = "none"),
##'                         covtype = "Matern5_2", orthogonal = TRUE, f.sim = f.sim, df.sim = df.sim)
##'
##'print(cpara.Het.OGP <- model$cpara)
##'Info.mx <- computeInfo(model)
##'# standard error of calibration parameter estimator
##'sd.cpara <- sqrt(diag(Info.mx)[1])
##'print(LCL <- qnorm(0.025, cpara.Het.OGP, sd.cpara)) # lower bound
##'print(UCL <- qnorm(0.975, cpara.Het.OGP, sd.cpara)) # upper bound
##' @export
computeInfo <- function(object){
  # number of unique locations
  n <- nrow(object$X0)
  N <- length(object$Z)
  # number of parameters
  m <- length(object$cpara) + length(object$theta) + length(object$nu_hat) +
    length(object$theta_g) + length(object$g) + length(object$nu_hat_var) + length(object$Delta)
  # precompute
  U <- matrix(0, ncol = length(object$mult), nrow = sum(object$mult))
  for(i in 1:ncol(U)) U[c(1,cumsum(object$mult)+1)[i]:cumsum(object$mult)[i],i] <- rep(1,object$mult[i])
  K <- cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype,
               orthogonal = object$orthogonal, f.sim = object$f.sim, df.sim = object$df.sim,
               cpara = object$cpara, MC.num = object$MC.num, inputBounds = object$inputBounds)
  nmean <- drop(rowSums(object$Kgi) %*% object$Delta / sum(object$Kgi))
  Kg <- cov_gen(X1 = object$X0, theta = object$theta_g, type = object$covtype)
  Lambda <- drop(nmean + Kg %*% (object$Kgi %*% (object$Delta - nmean)))
  Lambda <- exp(Lambda)
  NLambda <- rep(Lambda, times = object$mult)
  #V <- object$nu_hat * U %*% (K + diag(Lambda)) %*% t(U)
  V <- object$nu_hat * (U%*%K%*%t(U) + diag(NLambda))
  Vinv <- 1/object$nu_hat * (diag(1/NLambda) -
                              diag(1/NLambda) %*% U %*% solve(solve(K) + t(U) %*% diag(1/NLambda) %*% U) %*% t(U) %*% diag(1/NLambda))
  An <- diag(object$mult)
  Ainv <- diag(1/object$mult)
  Kgi <- object$Kgi
  Gg <- Kg + diag(object$eps + object$g/object$mult)
  W <- object$nu_hat_var * Gg
  Winv <- solve(W)

  M <- add_diag(Kgi * (-object$eps - object$g / object$mult), rep(1, n))
  rSKgi <- rowSums(Kgi)
  sKgi <- sum(Kgi)
  KgiD <- Kgi %*% (object$Delta - nmean)
  rsM <- rowSums(M)
  R <- diag(1,n) - matrix(1,n,n) %*% Kgi / sKgi


  ### compute derivatives of V with all parameters
  para_count <- 0
  dV <- array(0, dim = c(N, N, m))
  dW <- array(0, dim = c(n, n, m))
  ddW <- array(0, dim = c(n, n, m, m))

  # dV/dtheta
  for(i in 1:length(object$cpara)){
    dK_dtheta <- partial_cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype, arg = "cpara",
                                 orthogonal = object$orthogonal, f.sim = object$f.sim, df.sim = object$df.sim,
                                 cpara = object$cpara, MC.num = object$MC.num, inputBounds = object$inputBounds)[[i]]
    dV[,, i] <- object$nu_hat * U %*% dK_dtheta %*% t(U)
  }
  para_count <- para_count + length(object$cpara)
  # dV/dpsi
  if(length(object$theta) == 1){
    dK_dpsi <- partial_cov_gen(X1 = object$X0, theta = object$theta, arg = "theta_k", type = object$covtype,
                               orthogonal = object$orthogonal, f.sim = object$f.sim, df.sim = object$df.sim,
                               cpara = object$cpara, MC.num = object$MC.num, inputBounds = object$inputBounds)
    dV[,, para_count+1] <- object$nu_hat * U %*% dK_dpsi %*% t(U)
  }else{
    for(i in 1:length(object$theta)){
      dK_dpsi <- partial_cov_gen(X1 = object$X0[, i, drop = FALSE], theta = object$theta[i], arg = "theta_k", type = object$covtype,
                                 orthogonal = object$orthogonal, f.sim = object$f.sim, df.sim = object$df.sim,
                                 cpara = object$cpara, MC.num = object$MC.num, inputBounds = object$inputBounds)
      dV[,, para_count+i] <- object$nu_hat * U %*% dK_dpsi %*% t(U)
    }
  }
  para_count <- para_count + length(object$theta)
  #dV/dnu
  dV[,, para_count+1] <- V/object$nu_hat
  para_count <- para_count + 1
  #dV/dphi, dW/dphi, ddW/dphidphi
  if(length(object$theta_g) == 1){
    dKg_dphi <- partial_cov_gen(X1 = object$X0, theta = object$theta_g, arg = "theta_k", type = object$covtype) * Kg
    dV[,, para_count+1] <- object$nu_hat * diag(rep(Lambda * drop(dKg_dphi %*% KgiD - M %*% (dKg_dphi %*% KgiD) -
                                                                    (1 - rsM) * drop(rSKgi %*% dKg_dphi %*% (Kgi %*% object$Delta) * sKgi - rSKgi %*% object$Delta * (rSKgi %*% dKg_dphi %*% rSKgi))/sKgi^2), times = object$mult))
    dW[,, para_count+1] <- object$nu_hat_var * dKg_dphi
  }else{
    for(i in 1:length(object$theta)){
      dKg_dphi <- partial_cov_gen(X1 = object$X0[, i, drop = FALSE], theta = object$theta_g[i], arg = "theta_k", type = object$covtype) * Kg
      dV[,, para_count+i] <- object$nu_hat * diag(rep(Lambda * drop(dKg_dphi %*% KgiD - M %*% (dKg_dphi %*% KgiD) -
                                                                      (1 - rsM) * drop(rSKgi %*% dKg_dphi %*% (Kgi %*% object$Delta) * sKgi - rSKgi %*% object$Delta * (rSKgi %*% dKg_dphi %*% rSKgi))/sKgi^2), times = object$mult))
      dW[,, para_count+i] <- object$nu_hat_var * dKg_dphi
    }
  }
  if(length(object$theta_g) == 1){
    ddKg_dphiphi <- ddKg_ddphi_func(object)
    ddW[,,para_count+1,para_count+1] <- object$nu_hat_var * ddKg_dphiphi
  }else{
    for(i in 1:length(object$theta_g)){
      for(j in i:length(object$theta_g)){
        if(i == j) {
          ddW[,,para_count+i,para_count+i] <- object$nu_hat_var * ddKg_ddphi_func(object, arg = i) * cov_gen(X1 = object$X0[,-i,drop=FALSE], theta = object$theta_g[-i], type = object$covtype)
        }else{
          dlogKg_dphi_i <- partial_cov_gen(X1 = object$X0[, i, drop = FALSE], theta = object$theta_g[i], arg = "theta_k", type = object$covtype)
          dlogKg_dphi_j <- partial_cov_gen(X1 = object$X0[, j, drop = FALSE], theta = object$theta_g[j], arg = "theta_k", type = object$covtype)
          ddW[,,para_count+i,para_count+j] <- ddW[,,para_count+j,para_count+i] <- object$nu_hat_var * dlogKg_dphi_i * dlogKg_dphi_j * Kg
        }
      }
    }
  }
  para_count <- para_count + length(object$theta_g)
  #dV/dg, dW/dg, ddW/dphidg = 0, ddW/dgdg = 0
  dV[,, para_count+1] <- - object$nu_hat * diag(rep(Lambda * drop(-M %*% (KgiD/object$mult) - (1 - rsM) * drop(object$Delta %*% (Kgi %*% (rSKgi/object$mult)) * sKgi - rSKgi %*% object$Delta * sum(rSKgi^2/object$mult))/sKgi^2), times = object$mult))
  dW[,, para_count+1] <- object$nu_hat_var * Ainv
  dGinv_di <- -Winv %*% dW[,,para_count+1] %*% Winv * object$nu_hat_var

  para_count <- para_count + 1
  #dV/dnug, dW/dnug, ddW/dphidnug, ddW/dgdnug, ddW/dnugdnug
  dV[,, para_count+1] <- matrix(0, N, N)
  dW[,, para_count+1] <- Gg
  if(length(object$theta_g) == 1){ #ddW/dphidnug
    dKg_dphi <- partial_cov_gen(X1 = object$X0, theta = object$theta_g, arg = "theta_k", type = object$covtype) * Kg
    ddW[,,para_count-1,para_count+1] <- ddW[,,para_count+1,para_count-1] <- dKg_dphi
  }else{
    for(i in 1:length(object$theta)){
      dKg_dphi <- partial_cov_gen(X1 = object$X0[, i, drop = FALSE], theta = object$theta_g[i], arg = "theta_k", type = object$covtype) * Kg
      ddW[,,(para_count-length(object$theta_g)-1):(para_count-1),para_count+1] <-
        ddW[,,para_count+1,(para_count-length(object$theta_g)-1):(para_count-1)] <- dKg_dphi
    }
  }
  ddW[,,para_count,para_count+1] <- ddW[,,para_count+1,para_count] <- Ainv #ddW/dgdnug

  para_count <- para_count + 1
  #dV/dDelta
  for(i in 1:n){
    ei <- rep(0, n)
    ei[i] <- 1
    dV[,, para_count+i] <- object$nu_hat * diag(rep(Lambda * drop(Kg %*% Kgi %*% ei), times = object$mult))
  }

  H.base <- matrix(0, m, m)
  for(i in 1:m){
    for(j in i:m){
      H.base[i,j] <- H.base[j,i] <- -sum(diag(Vinv %*% dV[,,i] %*% Vinv %*% dV[,,j]))/2
    }
  }

  H.add <- matrix(0, m, m)
  # theta
  m1 <- length(object$cpara)
  if(length(object$cpara) == 1){
    H.add[1,1] <- drop(-object$df.sim(U%*%object$X0, object$cpara) %*% Vinv %*% object$df.sim(U%*%object$X0, object$cpara))
  }else{
    for(i in 1:length(object$cpara)){
      for(j in i:length(object$cpara)){
        H.add[i,j] <- H.add[j,i] <- drop(- object$df.sim(U%*%object$X0, object$cpara)[,i] %*% Vinv %*% object$df.sim(U%*%object$X0, object$cpara)[,j])
      }
    }
  }

  # psi and nu
  m2 <- length(c(object$theta, object$nu_hat))

  # phi, g, nug
  m3 <- length(c(object$theta_g, object$g, object$nu_hat_var))
  for(i in (m1+m2+1):(m1+m2+m3)){
    for(j in (m1+m2+1):(m1+m2+m3)){
      dWinv_di <- -Winv %*% dW[,,i] %*% Winv
      dWinv_dj <- -Winv %*% dW[,,j] %*% Winv
      ddWinv_didj <- Winv%*%(dW[,,i]%*%Winv%*%dW[,,j]+dW[,,j]%*%Winv%*%dW[,,i]-ddW[,,i,j])%*%Winv

      dGinv_di <- object$nu_hat_var * dWinv_di
      dGinv_dj <- object$nu_hat_var * dWinv_dj
      ddGinv_didj <- object$nu_hat_var * ddWinv_didj
      dR_i <-  sum(dGinv_di)*(matrix(1,n,n) %*% Kgi)/sKgi^2 - (matrix(1,n,n) %*% dGinv_di)/sKgi
      dR_j <-  sum(dGinv_dj)*(matrix(1,n,n) %*% Kgi)/sKgi^2 - (matrix(1,n,n) %*% dGinv_dj)/sKgi
      ddR_ij <- -2 * sum(dGinv_di) * sum(dGinv_dj)/sKgi^3 * (matrix(1,n,n) %*% Kgi) +
        1/sKgi^2 *(sum(ddGinv_didj)*(matrix(1,n,n) %*% Kgi) + sum(dGinv_dj)*(matrix(1,n,n) %*% dGinv_di) + sum(dGinv_di)*(matrix(1,n,n) %*% dGinv_dj)) -
        1/sKgi*(matrix(1,n,n) %*% ddGinv_didj)

      A1 <- t(R) %*% ddWinv_didj %*% R + t(dR_i) %*% dWinv_dj %*% R + t(R) %*% dWinv_dj %*% dR_i
      A2 <- t(dR_j) %*% dWinv_di %*% R + t(ddR_ij) %*% Winv %*% R + t(dR_j) %*% Winv %*% dR_i
      A3 <- t(A2)
      H.add[i,j] <- H.add[j,i] <-
        - sum(diag(Winv %*% ddW[,,i,j] + dWinv_di %*% dW[,,j]))/2 - (object$Delta - nmean) %*% (A1+A2+A3) %*% (object$Delta - nmean)/2
    }
  }
  # Delta
  for(j in (m1+m2+1):(m1+m2+m3)){
    dWinv_dj <- -Winv %*% dW[,,j] %*% Winv
    dGinv_dj <- object$nu_hat_var * dWinv_dj
    dR_j <-  sum(dGinv_dj)*(matrix(1,n,n) %*% Kgi)/sKgi^2 - (matrix(1,n,n) %*% dGinv_dj)/sKgi

    H.add[(m1+m2+m3+1):(m1+m2+m3+n),j] <- H.add[j,(m1+m2+m3+1):(m1+m2+m3+n)] <-
      - (t(R)%*% dWinv_dj %*% R + t(dR_j)%*% Winv %*% R + t(R) %*% Winv %*% dR_j) %*% object$Delta
  }

  for(i in 1:n){
    ei <- rep(0, n)
    ei[i] <- 1
    for(j in (m1+m2+1):(m1+m2+m3)){
      dWinv_dj <- -Winv %*% dW[,,j] %*% Winv
      H.add[m1+m2+m3+i,j] <- H.add[j,m1+m2+m3+i] <- - drop(ei %*% dWinv_dj %*% (object$Delta - nmean))
    }
  }
  H.add[(m1+m2+m3+1):(m1+m2+m3+n),(m1+m2+m3+1):(m1+m2+m3+n)] <- - Winv + outer(rSKgi, rSKgi)/sKgi/object$nu_hat_var

  H <- H.base + H.add
  colnames(H) <- c(ifelse(length(object$cpara) == 1, "cpara", paste0("cpara_",1:length(object$cpara))),
                   ifelse(length(object$theta) == 1, "theta",paste0("theta_",1:length(object$theta))),
                   "nu_hat",
                   ifelse(length(object$theta_g) == 1, "theta_g", paste0("theta_g_",1:length(object$theta_g))),
                   "g", "nu_hat_var",
                   paste0("Delta_",1:length(object$Delta)))
  H <- H[-c(m1+m2, m1+m2+m3),-c(m1+m2, m1+m2+m3)]
  return(solve(-H))
}
