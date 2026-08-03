
# KALMAN FILTER #

CTCRW_filter <- function(y, Hmat, beta1_vec, beta2_vec,
                         s1, s2, s3, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT_R(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ_R(beta1_vec[i], beta2_vec[i], s1, s2, s3, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    Z_i <- Z_full
    H_i <- diag(Hmat[i,], 3)
    
    v <- y[i,] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,]  <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}






###############################

# (using s_horiz and s_vert)

# DELETE
CTCRW_filter10 <- function(y, Hmat, beta1_vec, beta2_vec,
                           s_horiz, s_vert, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ(beta1_vec[i], beta2_vec[i], s_horiz, s_vert, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    Z_i <- Z_full
    H_i <- diag(Hmat[i,], 3)
    
    v <- y[i,] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,]  <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}
















# USE THIS VERSION AND CALL IT "CTCRW_filter1" for noerror, constanterror, and linearhdop
# simulation code files.
# change name back to CTCRW_filter1
CTCRW_filter111 <- function(y, Hmat, beta1_vec, beta2_vec,
                            s_horiz, s_vert, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ(beta1_vec[i], beta2_vec[i], s_horiz, s_vert, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    # mask missing observations
    obs_mask <- !is.na(y[i, ])
    
    # if all three are missing, just propagate without update
    if (!any(obs_mask)) {
      aest <- a_pred
      Pest <- P_pred
      a_f[i, ]  <- as.numeric(aest)
      P_f[[i]]  <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask, , drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- y[i, obs_mask] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i, ]  <- as.numeric(aest)
    P_f[[i]]  <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}







########################################################




# used in SpermWhale_LinearHDOP as "CTCRW_filter1

CTCRW_filter11 <- function(y, Hmat, beta1_vec, beta2_vec,
                           s_horiz, s_vert, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ(beta1_vec[i], beta2_vec[i], s_horiz, s_vert, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    # mask missing observations
    obs_mask <- !is.na(y[i, ])
    
    # if all three are missing, skip update
    if (!any(obs_mask)) {
      aest <- a_pred
      Pest <- P_pred
      a_f[i, ]  <- as.numeric(aest)
      P_f[[i]]  <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask, , drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- y[i, obs_mask] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i, ]  <- as.numeric(aest)
    P_f[[i]]  <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}




















CTCRW_filter1 <- function(y, Hmat, beta1_vec, beta2_vec,
                          s_horiz, s_vert, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ(beta1_vec[i], beta2_vec[i], s_horiz, s_vert, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    # mask missing observations
    obs_mask <- !is.na(y[i, ])
    
    # if all three are missing, skip update
    if (!any(obs_mask)) {
      aest <- a_pred
      Pest <- P_pred
      a_f[i, ]  <- as.numeric(aest)
      P_f[[i]]  <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask, , drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- y[i, obs_mask] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i, ]  <- as.numeric(aest)
    P_f[[i]]  <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}





#########################################







neg_loglikelihood_LinearError2 <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])   # horizontal process noise SD
  sigma2 <- exp(params["sigma2"])   # vertical process noise SD
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  delta <- rep(dt, nrow(data_aug))
  
  # 2-sigma process noise
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  # Build H matrix (linear HDOP)
  var0_xy  <- 0
  var1_xy  <- 0.1
  sd_depth <- 10
  
  Hmat <- build_Hmat_LinearError2(data_aug, var0_xy, var1_xy, sd_depth)
  
  # Initial state
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  P <- diag(6) * 1e2
  
  # 2-sigma filter
  filt <- CTCRW_filter1(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz   = s_horiz,
    s_vert    = s_vert,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  -filt$ll
}
















