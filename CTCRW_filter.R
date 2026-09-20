


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

























