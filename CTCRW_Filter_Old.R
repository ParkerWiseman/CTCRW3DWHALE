


# good code used before September 16, 2026`



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











###########################################



# *good code used before Sept 17, 2026


# KALMAN FILTER #




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














###############

# September 23, 2026




############################################################
# CTCRW_filter.R
#
# Kalman filter for the 3D CTCRW model
#
# State:
#   [x, vx, y, vy, depth, vdepth]
#
# The first observation is treated as being at the initial
# state time. Therefore:
#
#   delta[1] = 0
#   delta[i] = time[i] - time[i-1], i >= 2
#
# Missing x/y/depth observations are handled individually.
############################################################

CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  ##########################################################
  # Observation matrix
  #
  # observations:
  #   x     -> state 1
  #   y     -> state 3
  #   depth -> state 5
  ##########################################################
  
  Z_full <- matrix(0, 3, 6)
  
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  ##########################################################
  # Storage
  ##########################################################
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  ##########################################################
  # Initial state
  ##########################################################
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  ##########################################################
  # Kalman filter
  ##########################################################
  
  for (i in seq_len(N)) {
    
    ########################################################
    # PREDICTION
    ########################################################
    
    if (i == 1) {
      
      # Treating the first observation as the initial state,
      # instead of treating it as a movement step
      
      a_pred <- aest
      P_pred <- Pest
      
    } else {
      
      Tmat <- makeT(
        b1 = beta1_vec[i],
        b2 = beta2_vec[i],
        delta = delta[i]
      )
      
      Qmat <- makeQ(
        b1 = beta1_vec[i],
        b2 = beta2_vec[i],
        s_horiz = s_horiz,
        s_vert = s_vert,
        delta = delta[i]
      )
      
      a_pred <- as.numeric(Tmat %*% aest)
      
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    }
    
    ########################################################
    # (major stability improvement)
    ########################################################
    
    P_pred <- (P_pred + t(P_pred)) / 2
    
    ########################################################
    # Save prediction
    ########################################################
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    ########################################################
    # Determine which observations are available
    ########################################################
    
    obs_mask <- !is.na(y[i,])
    
    ########################################################
    # If everything is missing, simply propagate the state
    ########################################################
    
    if (!any(obs_mask)) {
      
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      
      next
    }
    
    ########################################################
    # Observation matrix for available observations
    ########################################################
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    
    H_i <- diag(Hmat[i, obs_mask],sum(obs_mask))
    
    ########################################################
    # Innovation
    ########################################################
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    ########################################################
    # Innovation covariance
    ########################################################
    # Old Code: Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    
    Fmat <- (Fmat + t(Fmat)) / 2
    
    ########################################################
    # Using Cholesky decomposition to compute inverse of F matrix,
    # which is more numerically stable than directly using solve() and det().
    # Old Code: invF <- solve(Fmat)
    ########################################################
    
    cholF <- tryCatch(chol(Fmat), error = function(e) NULL)
    
    ########################################################
    # If covariance matrix is not positive definite, then we return a failed likelihood
    # (i.e., it won't just crash if Fmat is not invertible) 
    ########################################################
    
    if (is.null(cholF)) {
      
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    ########################################################
    # Log determinant
    ########################################################
    
    logdetF <- 2 * sum(log(diag(cholF)))
    
    ########################################################
    # Quadratic likelihood term
    ########################################################
    
    F_inv_v <- backsolve(cholF, forwardsolve(cholF, v))
    
    ll <- ll - 0.5 * (logdetF + sum(v * F_inv_v))
    
    ########################################################
    # Kalman gain
    #
    # Old filter used: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    # Here I use: K <- P_pred %*% t(Z_i) %*% F_inv
    ########################################################
    
    F_inv <- chol2inv(cholF)
    
    K <- P_pred %*% t(Z_i) %*% F_inv
    
    ########################################################
    # State update
    ########################################################
    
    aest <- a_pred + K %*% v
    
    ########################################################
    # Old code used: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    # This code uses the Joseph-form covariance update:
    # Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    # This form of update guarantees: symmetry, positive‑definiteness, numerical stability
    ########################################################
    
    I6 <- diag(6)
    
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    Pest <- (Pest + t(Pest)) / 2
    
    ########################################################
    # Save filtered state
    ########################################################
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  ##########################################################
  # Return
  ##########################################################
  
  list(
    ll = ll,
    a_f = a_f,
    P_f = P_f,
    a_p = a_p,
    P_p = P_p
  )
}















##############

# evening of September 23


CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    Fmat <- (Fmat + t(Fmat)) / 2
    
    ########################################################
    # Using solve() instead of Cholesky
    ########################################################
    
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # F_inv <- invF
    
    #K <- P_pred %*% t(Z_i) %*% invF
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    #I6 <- diag(6)
    #
    #Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) +
    #  K %*% H_i %*% t(K)
    #
    #Pest <- (Pest + t(Pest)) / 2
    
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}











########################################

# night of September 23

#> cat("Horizontal autocorrelation time tau1 (seconds): ", tau1_hat, "\n")
#Horizontal autocorrelation time tau1 (seconds):  7062.616 
#> cat("Horizontal autocorrelation time tau1 (hours):   ", tau1_hat / 3600, "\n")
#Horizontal autocorrelation time tau1 (hours):    1.961838 
#> cat("Vertical autocorrelation time tau2 (seconds):   ", tau2_hat, "\n")
#Vertical autocorrelation time tau2 (seconds):    171.2704 
#> cat("Vertical autocorrelation time tau2 (hours):     ", tau2_hat / 3600, "\n")
#Vertical autocorrelation time tau2 (hours):      0.04757511 
#> cat("Horizontal RMS speed nu1 (m/s):                 ", nu1_hat, "\n")
#Horizontal RMS speed nu1 (m/s):                  0.753836 
#> cat("Vertical RMS speed nu2 (m/s):                   ", nu2_hat, "\n")
#Vertical RMS speed nu2 (m/s):                    0.4362287 








CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  Tmat <- diag(6)
  Qmat <- matrix(0,6,6)
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    Fmat <- (Fmat + t(Fmat)) / 2
    
    ########################################################
    # Using solve() instead of Cholesky
    ########################################################
    
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # F_inv <- invF
    
    #K <- P_pred %*% t(Z_i) %*% invF
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    #I6 <- diag(6)
    #
    #Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) +
    #  K %*% H_i %*% t(K)
    #
    #Pest <- (Pest + t(Pest)) / 2
    
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}








######################################################
######################################################
######################################################


# Old Code
# Night of September 28, 2026



CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      Tmat <- diag(6)
      Qmat <- matrix(0, 6, 6)
      a_pred <- aest
      P_pred <- Pest
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    ########################################################
    # Using solve() instead of Cholesky
    ########################################################
    
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # F_inv <- invF
    
    ######K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    ######Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Pest <- (Pest + t(Pest)) / 2
    
    #Old filter used this Pest update: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}

































CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    ########################################################
    # Using solve() instead of Cholesky
    ########################################################
    
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # F_inv <- invF
    
    K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Pest <- (Pest + t(Pest)) / 2
    
    #Old filter used this Pest update: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}














###############################



CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    ########################################################
    # Using solve() instead of Cholesky
    ########################################################
    
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # F_inv <- invF
    
    K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Pest <- (Pest + t(Pest)) / 2
    
    #Old filter used this Pest update: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}














#####################################################




### Newest code that still uses the RTS smoother ###



############################################################
# CTCRW_filter.R
############################################################

CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  # NEW: store Tmat for smoother
  T_array <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
      
      # First Tmat is identity (no movement at first observation)
      Tmat <- diag(6)
      T_array[[i]] <- Tmat
      
    } else {
      Tmat <- makeT(
        b1 = beta1_vec[i],
        b2 = beta2_vec[i],
        delta = delta[i]
      )
      
      Qmat <- makeQ(
        b1 = beta1_vec[i],
        b2 = beta2_vec[i],
        s_horiz = s_horiz,
        s_vert = s_vert,
        delta = delta[i]
      )
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
      
      # Store Tmat
      T_array[[i]] <- Tmat
    }
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # Standard Kalman gain
    K <- P_pred %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    # Joseph-form covariance update
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) +
      K %*% H_i %*% t(K)
    
    Pest <- (Pest + t(Pest)) / 2
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(
    ll = ll,
    a_f = a_f,
    P_f = P_f,
    a_p = a_p,
    P_p = P_p,
    T_array = T_array   # smoother will use this
  )
}



















###########################################################################



############################################################
# CTCRW_filter.R — 3D CTCRW, D&K-compatible
############################################################

CTCRW_filter1 <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  # 3D observation: x, y, depth
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  # Store for D&K smoother
  T_array <- vector("list", N)
  L_array <- vector("list", N)
  F_array <- vector("list", N)
  v_list  <- vector("list", N)
  Z_list  <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      Tmat <- diag(6)
      Qmat <- matrix(0, 6, 6)
      a_pred <- aest
      P_pred <- Pest
    } else {
      Tmat <- makeT(
        b1 = beta1_vec[i],
        b2 = beta2_vec[i],
        delta = delta[i]
      )
      
      Qmat <- makeQ(
        b1 = beta1_vec[i],
        b2 = beta2_vec[i],
        s_horiz = s_horiz,
        s_vert = s_vert,
        delta = delta[i]
      )
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    }
    
    T_array[[i]] <- Tmat
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      # No observations at this time
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      
      Z_list[[i]] <- Z_full[FALSE,,drop=FALSE]  # 0 x 6
      v_list[[i]] <- numeric(0)
      F_array[[i]] <- matrix(0, 0, 0)
      L_array[[i]] <- Tmat
      
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    
    Z_list[[i]] <- Z_i
    v_list[[i]] <- v
    F_array[[i]] <- Fmat
    
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    if (is.null(invF)) {
      return(list(ll = -Inf))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    K <- P_pred %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    L_array[[i]] <- Tmat - K %*% Z_i
    
    I6 <- diag(6)
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) +
      K %*% H_i %*% t(K)
    
    Pest <- (Pest + t(Pest)) / 2
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(
    ll      = ll,
    a_f     = a_f,
    P_f     = P_f,
    a_p     = a_p,
    P_p     = P_p,
    T_array = T_array,
    L_array = L_array,
    F_array = F_array,
    v_list  = v_list,
    Z_list  = Z_list
  )
}




###############################################




############################################################
# CTCRW_filter.R
############################################################

CTCRW_filter <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  # NEW: store Tmat for smoother instead of recalculating T during the smoothing step
  T_array <- vector("list", N)
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
      
      # First Tmat is identity (no movement at first observation)
      Tmat <- diag(6)
      T_array[[i]] <- Tmat
      
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
      
      # Store Tmat
      T_array[[i]] <- Tmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    # Using solve() instead of Cholesky
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    ##### Standard Kalman gain
    ###K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain:
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    ###I6 <- diag(6)
    
    #### "Joseph-form" covariance update
    ###Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Old filter used this Pest update:
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    # Pest <- (Pest + t(Pest)) / 2
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array)
}


















############################################################
# CTCRW_filter.R
############################################################

CTCRW_filter <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  # NEW: store Tmat for smoother instead of recalculating T during the smoothing step
  T_array <- vector("list", N)
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
      
      # First Tmat is identity (no movement at first observation)
      Tmat <- diag(6)
      T_array[[i]] <- Tmat
      
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
      
      # Store Tmat
      T_array[[i]] <- Tmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    # Using solve() instead of Cholesky
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # Standard Kalman gain
    K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    # "Joseph-form" covariance update
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Old filter used this Pest update: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    # Pest <- (Pest + t(Pest)) / 2
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array)
}














##############################################





############################################################
# CTCRW_filter.R
############################################################

CTCRW_filter <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  # NEW: store Tmat for smoother instead of recalculating T during the smoothing step
  T_array <- vector("list", N)
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
      
      # First Tmat is identity (no movement at first observation)
      Tmat <- diag(6)
      T_array[[i]] <- Tmat
      
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
      
      # Store Tmat
      T_array[[i]] <- Tmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    # Using solve() instead of Cholesky
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # Standard Kalman gain
    K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    # "Joseph-form" covariance update
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Old filter used this Pest update: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    # Pest <- (Pest + t(Pest)) / 2
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array)
}












#####################################################


# October 1, 2026



############################################################
# CTCRW_filter.R
############################################################

CTCRW_filter <- function(
    y,
    Hmat,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  # NEW: store Tmat for smoother instead of recalculating T during the smoothing step
  T_array <- vector("list", N)
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
      
      # First Tmat is identity (no movement at first observation)
      Tmat <- diag(6)
      T_array[[i]] <- Tmat
      
    } else {
      Tmat <- makeT(b1 = beta1_vec[i], b2 = beta2_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta1_vec[i], b2 = beta2_vec[i],
                    s_horiz = s_horiz, s_vert = s_vert, delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
      
      # Store Tmat
      T_array[[i]] <- Tmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    # Using solve() instead of Cholesky
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # Standard Kalman gain
    K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    # "Joseph-form" covariance update
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Old filter used this Pest update: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    # Pest <- (Pest + t(Pest)) / 2
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array)
}













####################################################################################


############################################################
# CTCRW_filter.R
############################################################

CTCRW_filter <- function(
    y,
    Hmat,
    beta_horiz_vec,
    beta_vert_vec,
    var_horiz_vec,
    var_vert_vec,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA_real_, N, 6)
  P_f <- vector("list", N)
  
  a_p <- matrix(NA_real_, N, 6)
  P_p <- vector("list", N)
  
  aest <- as.numeric(a)
  Pest <- P
  
  # NEW: store Tmat for smoother instead of recalculating T during the smoothing step
  T_array <- vector("list", N)
  
  ll <- 0
  
  for (i in seq_len(N)) {
    
    if (i == 1) {
      a_pred <- aest
      P_pred <- Pest
      
      # First Tmat is identity (no movement at first observation)
      Tmat <- diag(6)
      T_array[[i]] <- Tmat
      
    } else {
      Tmat <- makeT(b1 = beta_horiz_vec[i], b2 = beta_vert_vec[i], delta = delta[i])
      Qmat <- makeQ(b1 = beta_horiz_vec[i], b2 = beta_vert_vec[i],
                    var_horiz = var_horiz_vec[i], var_vert = var_vert_vec[i], delta = delta[i])
      
      a_pred <- as.numeric(Tmat %*% aest)
      P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
      
      # Store Tmat
      T_array[[i]] <- Tmat
    }
    
    #P_pred <- (P_pred + t(P_pred)) / 2
    
    a_p[i,] <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_mask <- !is.na(y[i,])
    
    if (!any(obs_mask)) {
      aest <- as.numeric(a_pred)
      Pest <- P_pred
      a_f[i,] <- aest
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_mask,,drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    v <- as.numeric(y[i, obs_mask] - Z_i %*% a_pred)
    
    Fmat <- Z_i %*% P_pred %*% t(Z_i) + H_i
    #Fmat <- (Fmat + t(Fmat)) / 2
    
    # Using solve() instead of Cholesky
    invF <- tryCatch(solve(Fmat), error = function(e) NULL)
    
    if (is.null(invF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    logdetF <- tryCatch(log(det(Fmat)), error = function(e) -Inf)
    
    if (!is.finite(logdetF)) {
      return(list(ll = -Inf, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array))
    }
    
    ll <- ll - 0.5 * (logdetF + t(v) %*% invF %*% v)
    
    # Standard Kalman gain
    K <- P_pred %*% t(Z_i) %*% invF
    
    #Old filter used this Kalman gain: K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    
    I6 <- diag(6)
    
    # "Joseph-form" covariance update
    Pest <- (I6 - K %*% Z_i) %*% P_pred %*% t(I6 - K %*% Z_i) + K %*% H_i %*% t(K)
    
    #Old filter used this Pest update: Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    # Pest <- (Pest + t(Pest)) / 2
    
    a_f[i,] <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p, T_array = T_array)
}










#######################################################################

# Morning of October 3, 2026



############################################################
# CTCRW_filter.R — 3D version of instructor’s filter
############################################################

CTCRW_filter <- function(
    y,
    Hmat,
    beta_horiz_vec,
    beta_vert_vec,
    var_horiz_vec,
    var_vert_vec,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  # Observation matrix (3x6)
  Z <- matrix(0, 3, 6)
  Z[1,1] <- 1   # x
  Z[2,3] <- 1   # y
  Z[3,5] <- 1   # depth
  
  # Storage
  aest <- matrix(0, N, 6)                 # filtered means
  Pest <- array(0, c(6, 6, N))    # filtered covariances
  
  a_pred <- matrix(0, N, 6)               # predicted means
  P_pred <- array(0, c(6, 6, N))  # predicted covariances
  
  v <- matrix(0, N, 3)                      # innovations (full 3D)
  F <- array(0, c(3, 3, N))           # innovation covariances (3x3)
  K <- array(0, c(6, 3, N))         # Kalman gains (6x3)
  L <- array(0, c(6, 6, N))       # L matrices (6x6)
  T_array <- array(0, c(6, 6, N)) # transition matrices (6x6)
  
  # Initialize
  aest[1,] <- a
  Pest[,,1] <- P
  
  ll <- 0
  
  # Loop over time
  for (i in 1:(N-1)) {
    
    # Transition matrix
    Tmat <- makeT(
      b1 = beta_horiz_vec[i+1],
      b2 = beta_vert_vec[i+1],
      delta = delta[i+1]
    )
    T_array[,,i] <- Tmat
    
    # Process noise
    Qmat <- makeQ(
      b1 = beta_horiz_vec[i+1],
      b2 = beta_vert_vec[i+1],
      var_horiz = var_horiz_vec[i+1],
      var_vert = var_vert_vec[i+1],
      delta = delta[i+1]
    )
    
    # Prediction
    a_pred[i+1,] <- Tmat %*% aest[i,]
    P_pred[,,i+1] <- Tmat %*% Pest[,,i] %*% t(Tmat) + Qmat
    
    # Missing observation?
    obs_mask <- !is.na(y[i+1,])
    
    if (!any(obs_mask)) {
      aest[i+1,] <- a_pred[i+1,]
      Pest[,,i+1] <- P_pred[,,i+1]
      L[,,i] <- Tmat
      next
    }
    
    # Subset Z and H to observed components
    Z_i <- Z[obs_mask,,drop=FALSE]
    H_i <- diag(Hmat[i+1, obs_mask], sum(obs_mask))
    
    # Innovation only for observed components
    v_i <- y[i+1, obs_mask] - Z_i %*% a_pred[i+1,]
    
    # Innovation covariance
    F_i <- Z_i %*% P_pred[,,i+1] %*% t(Z_i) + H_i
    
    invF_i <- solve(F_i)
    
    # Kalman gain for observed components
    K_i <- P_pred[,,i+1] %*% t(Z_i) %*% invF_i  # 6 x k
    
    # Embed v_i and K_i back into full 3D space
    v[i+1, obs_mask] <- v_i
    F[obs_mask, obs_mask, i+1] <- F_i
    
    # Full 6x3 K: zeros then fill observed columns
    K_full <- matrix(0, 6, 3)
    K_full[, obs_mask] <- K_i
    K[,,i+1] <- K_full
    
    # Update mean using full K and full v row
    aest[i+1,] <- a_pred[i+1,] + K_full %*% v[i+1,]
    
    # Update covariance (Joseph form) using full Z and full H
    H_full <- diag(Hmat[i+1,], 3)
    I6 <- diag(6)
    
    Pest[,,i+1] <- (I6 - K_full %*% Z) %*% P_pred[,,i+1] %*%
      t(I6 - K_full %*% Z) + K_full %*% H_full %*% t(K_full)
    
    # L matrix for smoother
    L[,,i] <- Tmat - K_full %*% Z
    
    # Log-likelihood (only observed components)
    ll <- ll - 0.5 * (log(det(F_i)) + t(v_i) %*% invF_i %*% v_i)
  }
  
  list(
    ll = ll,
    aest = aest,
    Pest = Pest,
    a_pred = a_pred,
    P_pred = P_pred,
    v = v,
    F = F,
    K = K,
    L = L,
    T_array = T_array
  )
}








############################################################################



############################################################
# CTCRW_filter.R — 3D version of instructor’s filter
############################################################

CTCRW_filter <- function(
    y,
    Hmat,
    beta_horiz_vec,
    beta_vert_vec,
    var_horiz_vec,
    var_vert_vec,
    delta,
    a,
    P
) {
  
  N <- nrow(y)
  
  # Observation matrix (3x6)
  Z <- matrix(0, 3, 6)
  Z[1,1] <- 1   # x
  Z[2,3] <- 1   # y
  Z[3,5] <- 1   # depth
  
  # Storage
  aest <- matrix(0, N, 6)                 # filtered means
  Pest <- array(0, c(6, 6, N))    # filtered covariances
  
  a_pred <- matrix(0, N, 6)               # predicted means
  P_pred <- array(0, c(6, 6, N))  # predicted covariances
  
  v <- matrix(0, N, 3)                      # innovations (full 3D)
  F <- array(0, c(3, 3, N))           # innovation covariances (3x3)
  K <- array(0, c(6, 3, N))         # Kalman gains (6x3)
  L <- array(0, c(6, 6, N))       # L matrices (6x6)
  T_array <- array(0, c(6, 6, N)) # transition matrices (6x6)
  
  # Initialize
  aest[1,] <- a
  Pest[,,1] <- P
  
  ll <- 0
  
  # Loop over time
  for (i in 1:(N-1)) {
    
    # Transition matrix
    Tmat <- makeT(
      b1 = beta_horiz_vec[i],
      b2 = beta_vert_vec[i],
      delta = delta[i+1]
    )
    T_array[,,i] <- Tmat
    
    
    
    # Process noise
    Qmat <- makeQ(
      b1 = beta_horiz_vec[i],
      b2 = beta_vert_vec[i],
      var_horiz = var_horiz_vec[i],
      var_vert = var_vert_vec[i],
      delta = delta[i+1]
    )
    
    
    
    
    # Prediction
    a_pred[i+1,] <- Tmat %*% aest[i,]
    P_pred[,,i+1] <- Tmat %*% Pest[,,i] %*% t(Tmat) + Qmat
    
    # Missing observation?
    obs_mask <- !is.na(y[i+1,])
    
    #if (!any(obs_mask)) {
    #  aest[i+1,] <- a_pred[i+1,]
    #  Pest[,,i+1] <- P_pred[,,i+1]
    #  L[,,i] <- Tmat
    #  next
    #}
    
    # Subset Z and H to observed components
    Z_i <- Z[obs_mask,,drop=FALSE]
    H_i <- diag(Hmat[i+1, obs_mask], sum(obs_mask))
    
    # Innovation only for observed components
    v_i <- y[i+1, obs_mask] - Z_i %*% a_pred[i+1,]
    
    # Innovation covariance
    F_i <- Z_i %*% P_pred[,,i+1] %*% t(Z_i) + H_i
    
    invF_i <- solve(F_i)
    
    # Kalman gain for observed components
    K_i <- P_pred[,,i+1] %*% t(Z_i) %*% invF_i  # 6 x k
    
    # Embed v_i and K_i back into full 3D space
    v[i+1, obs_mask] <- v_i
    F[obs_mask, obs_mask, i+1] <- F_i
    
    # Full 6x3 K: zeros then fill observed columns
    K_full <- matrix(0, 6, 3)
    K_full[, obs_mask] <- K_i
    K[,,i+1] <- K_full
    
    # Update mean using full K and full v row
    aest[i+1,] <- a_pred[i+1,] + K_full %*% v[i+1,]
    
    # Update covariance (Joseph form) using full Z and full H
    H_full <- diag(Hmat[i+1,], 3)
    I6 <- diag(6)
    
    Pest[,,i+1] <- (I6 - K_full %*% Z) %*% P_pred[,,i+1] %*%
      t(I6 - K_full %*% Z) + K_full %*% H_full %*% t(K_full)
    
    # L matrix for smoother
    L[,,i] <- Tmat - K_full %*% Z
    
    # Log-likelihood (only observed components)
    ll <- ll - 0.5 * (log(det(F_i)) + t(v_i) %*% invF_i %*% v_i)
  }
  
  list(
    ll = ll,
    aest = aest,
    Pest = Pest,
    a_pred = a_pred,
    P_pred = P_pred,
    v = v,
    F = F,
    K = K,
    L = L,
    T_array = T_array
  )
}






















