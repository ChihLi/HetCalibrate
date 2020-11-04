cov_ortho_gen <- function(X1, X2 = NULL, theta, type = c("Gaussian", "Matern5_2", "Matern3_2"), 
                          f.sim, df.sim = NULL, cpara, MC.num = 1000, inputBounds = NULL, nugget = 1e-6){
  
  if(is.null(inputBounds)){
    inputBounds <- find_reps(X1, rep(1, nrow(X1)), rescale = TRUE)$inputBounds
  }
  
  cov.ori <- eval(call(paste0("cov_",type), X1, X2, theta))
  
  #xi <- matrix(runif(MC.num, 0, 1), ncol = ncol(X1))
  if(ncol(X1) == 1){
    xi <- matrix(seq(0, 1, length.out = MC.num), ncol = ncol(X1))
  }else{
    xi <- sobol(MC.num, ncol(X1))
  }

  xi <- xi * matrix(rep(inputBounds[2,] - inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) + 
    matrix(rep(inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) 
  
  Wm <- eval(call(paste0("cov_",type), X1 = xi, X2 = NULL, theta))
  #wx1 <- apply(X1, 1, function(x) eval(call(paste0("cov_",type), as.matrix(x, ncol = ncol(X1)), xi, theta)))
  wx1 <- apply(X1, 1, function(x) eval(call(paste0("cov_",type), matrix(x, ncol = ncol(X1)), xi, theta)))
  if(is.null(X2)) {
    wx2 <- wx1
  }else{
    #wx2 <- apply(X2, 1, function(x) eval(call(paste0("cov_",type), as.matrix(x, ncol = ncol(X2)), xi, theta)))
    wx2 <- apply(X2, 1, function(x) eval(call(paste0("cov_",type), matrix(x, ncol = ncol(X2)), xi, theta)))
  }
  
  # if(is.null(df.sim)){
  #   Fm <- apply(xi, 1, function(xi.i) {
  #     f.tmp <- function(x) c(f.sim(xi.i, x))
  #     grad(f.tmp, cpara)
  #     })
  #   if(ncol(X1) == 1) Fm <- matrix(Fm, ncol = 1)
  # }else{
  Fm <- df.sim(xi, cpara)
  if(length(cpara) == 1) Fm <- matrix(Fm, ncol = 1)
  # }
  
  rm.index <- apply(Fm, 2, function(x) all(x==0))
  if(all(rm.index)){
    cov.ortho <- cov.ori
  }else{
    if(length(cpara) == 1){
      FWF_inverse <- 1/(t(Fm) %*% Wm %*% Fm)
    }else{
      if(any(rm.index)) {
        Fm <- Fm[,!rm.index]
      }
      FWF <- t(Fm) %*% Wm %*% Fm
      FWF_inverse <- solve(FWF + diag(nugget * diag(FWF)))
    }
    cov.ortho <- cov.ori - t(wx1) %*% Fm %*% FWF_inverse %*% t(Fm) %*% wx2
  }
  #FWF_inverse <- solve(t(Fm) %*% Wm %*% Fm + diag(delta, length(cpara)))

  
  #cov.ortho <- (cov.ortho + t(cov.ortho))/2 # ensure symmetric 
  
  return(cov.ortho)
}

partial_cov_ortho_gen <- function(X1, X2 = NULL, theta, type = c("Gaussian", "Matern5_2", "Matern3_2"), arg = "theta_k", 
                                  f.sim, df.sim = NULL, cpara, MC.num = 1000, inputBounds = NULL, nugget = 1e-6){
  
  if(is.null(inputBounds)){
    inputBounds <- find_reps(X1, rep(1, nrow(X1)), rescale = TRUE)$inputBounds
  }
  
  #xi <- matrix(runif(MC.num, 0, 1), ncol = ncol(X1))
  if(ncol(X1) == 1){
    xi <- matrix(seq(0, 1, length.out = MC.num), ncol = ncol(X1))
  }else{
    xi <- sobol(MC.num, ncol(X1))
  }
  xi <- xi * matrix(rep(inputBounds[2,] - inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) + 
    matrix(rep(inputBounds[1,], each = nrow(xi)), ncol = ncol(xi)) 
  
  Wm <- eval(call(paste0("cov_",type), X1 = xi, X2 = NULL, theta))
  #wx1 <- apply(X1, 1, function(x) eval(call(paste0("cov_",type), as.matrix(x, ncol = ncol(X1)), xi, theta)))
  wx1 <- apply(X1, 1, function(x) eval(call(paste0("cov_",type), matrix(x, ncol = ncol(X1)), xi, theta)))
  if(is.null(X2)) {
    wx2 <- wx1
  }else{
    #wx2 <- apply(X2, 1, function(x) eval(call(paste0("cov_",type), as.matrix(x, ncol = ncol(X2)), xi, theta)))
    wx2 <- apply(X2, 1, function(x) eval(call(paste0("cov_",type), matrix(x, ncol = ncol(X2)), xi, theta)))
  }
  # if(is.null(df.sim)){
  #   Fm <- apply(xi, 1, function(xi.i) {
  #     f.tmp <- function(x) c(f.sim(xi.i, x))
  #     grad(f.tmp, cpara)
  #   })
  #   if(ncol(X1) == 1) Fm <- matrix(Fm, ncol = 1)
  # }else{
  Fm <- df.sim(xi, cpara)
  if(length(cpara) == 1) Fm <- matrix(Fm, ncol = 1)
  # }
  
  if(arg == "theta_k"){
    #partial_wx1 <- apply(X1, 1, function(x) eval(call(paste0("partial_d_k_", type, "_dtheta_k"), as.matrix(x, ncol = ncol(X1)), xi, theta))) * wx1
    partial_wx1 <- apply(X1, 1, function(x) eval(call(paste0("partial_d_k_", type, "_dtheta_k"), matrix(x, ncol = ncol(X1)), xi, theta))) * wx1
    partial_Wm <- eval(call(paste0("partial_d_C_", type, "_dtheta_k"), X1 = xi, theta)) * Wm
    
    if(is.null(X2)) {
      partial_cov.ori <- eval(call(paste0("partial_d_C_", type, "_dtheta_k"), X1, theta)) * eval(call(paste0("cov_", type), X1, X2 = NULL, theta))
      partial_wx2 <- partial_wx1
    }else{
      partial_cov.ori <- eval(call(paste0("partial_d_k_", type, "_dtheta_k"), X1, X2, theta)) * eval(call(paste0("cov_", type), X1, X2, theta))
      #partial_wx2 <- apply(X2, 1, function(x) eval(call(paste0("partial_d_k_", type, "_dtheta_k"), as.matrix(x, ncol = ncol(X2)), xi, theta))) * wx2
      partial_wx2 <- apply(X2, 1, function(x) eval(call(paste0("partial_d_k_", type, "_dtheta_k"), matrix(x, ncol = ncol(X2)), xi, theta))) * wx2
    }
    
    rm.index <- apply(Fm, 2, function(x) all(x==0))
    if(all(rm.index)){
      partial_cov.ortho <- matrix(0, MC.num, MC.num)
    }else{
      if(length(cpara) == 1){
        FWF_inverse <- 1/(t(Fm) %*% Wm %*% Fm)
      }else{
        if(any(rm.index)) {
          Fm <- Fm[,!rm.index]
        }
        FWF <- t(Fm) %*% Wm %*% Fm
        FWF_inverse <- solve(FWF + diag(nugget * diag(FWF)))
      }
      partial_cov.ortho <- partial_cov.ori - 
        t(partial_wx1) %*% Fm %*% FWF_inverse %*% t(Fm) %*% wx2 - 
        t(wx1) %*% Fm %*% FWF_inverse %*% t(Fm) %*% partial_wx2 + 
        t(wx1) %*% Fm %*% FWF_inverse %*% (t(Fm) %*% partial_Wm %*% Fm) %*% FWF_inverse %*% t(Fm) %*% wx2
    }
    

  }else if(arg == "cpara"){
    partial_cov.ortho <- vector("list", length(cpara))
    for(i in 1:length(cpara)){
      ddf.sim <- function(x, cpara){
        # if(is.vector(x)) x <- matrix(x, ncol = 1)
        # apply(x, 1, function(y) {
        #   f.tmp <- function(para) c(df.sim(y, para))
        #   grad(f.tmp, cpara)
        # })
        if(length(cpara) == 1){
          f.tmp <- function(cpara.tmp, x.tmp) c(df.sim(x.tmp, cpara.tmp))
          return(gradient(f.tmp, cpara, x.tmp = x, centered = TRUE, pert = 1e-4))
        }else{
          dFm <- matrix(0, ncol = length(cpara), nrow = nrow(xi))
          for(j in 1:length(cpara)){
            f.tmp <- function(cpara.tmp, x.tmp) {
              tmp.vt <- cpara
              tmp.vt[i] <- cpara.tmp
              df.sim(x.tmp, tmp.vt)[,j]
            }
            dFm[,j] <- gradient(f.tmp, cpara[i], x.tmp = x, centered = TRUE, pert = 1e-4)
          }
          return(dFm)
        }
      }
      partial_Fm <- ddf.sim(xi, cpara)  
      if(length(cpara) == 1) partial_Fm <- matrix(partial_Fm, ncol = 1)
      
      
      rm.index <- apply(Fm, 2, function(x) all(x==0))
      if(all(rm.index)){
        partial_cov.ortho[[i]] <- matrix(0, MC.num, MC.num)
      }else{
        if(length(cpara) == 1){
          FWF_inverse <- 1/(t(Fm) %*% Wm %*% Fm)
          pFWF <- t(partial_Fm) %*% Wm %*% Fm
          partial_cov.ortho[[i]] <- t(wx1) %*% Fm %*% FWF_inverse %*% pFWF %*% FWF_inverse %*% t(Fm) %*% wx2 +
            t(wx1) %*% Fm %*% FWF_inverse %*% pFWF %*% FWF_inverse %*% t(Fm) %*% wx2 -
            t(wx1) %*% partial_Fm %*% FWF_inverse %*% t(Fm) %*% wx2 - t(wx1) %*% Fm %*% FWF_inverse %*% t(partial_Fm) %*% wx2
        }else{
          if(any(rm.index)) {
            Fm.tmp <- Fm[,!rm.index]
            partial_Fm.tmp <- partial_Fm[,!rm.index]
            FWF <- t(Fm.tmp) %*% Wm %*% Fm.tmp
            pFWF <- t(partial_Fm.tmp) %*% Wm %*% Fm.tmp
            FWF_inverse <- solve(FWF + diag(nugget * diag(FWF)))
            partial_cov.ortho[[i]] <- t(wx1) %*% Fm.tmp %*% FWF_inverse %*% pFWF %*% FWF_inverse %*% t(Fm.tmp) %*% wx2 +
              t(wx1) %*% Fm.tmp %*% FWF_inverse %*% pFWF %*% FWF_inverse %*% t(Fm.tmp) %*% wx2 -
              t(wx1) %*% partial_Fm.tmp %*% FWF_inverse %*% t(Fm.tmp) %*% wx2 - t(wx1) %*% Fm.tmp %*% FWF_inverse %*% t(partial_Fm.tmp) %*% wx2
          }else{
            FWF <- t(Fm) %*% Wm %*% Fm
            pFWF <- t(partial_Fm) %*% Wm %*% Fm
            FWF_inverse <- solve(FWF + diag(nugget * diag(FWF)))
            partial_cov.ortho[[i]] <- t(wx1) %*% Fm %*% FWF_inverse %*% pFWF %*% FWF_inverse %*% t(Fm) %*% wx2 +
              t(wx1) %*% Fm %*% FWF_inverse %*% pFWF %*% FWF_inverse %*% t(Fm) %*% wx2 -
              t(wx1) %*% partial_Fm %*% FWF_inverse %*% t(Fm) %*% wx2 - t(wx1) %*% Fm %*% FWF_inverse %*% t(partial_Fm) %*% wx2
          }
        }

      }
    }

    
  }
  
  
  #partial_cov.ortho <- (partial_cov.ortho + t(partial_cov.ortho))/2 # ensure symmetric 
  
  return(partial_cov.ortho)
}
