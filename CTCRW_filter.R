


############################################################
# CTCRW_filter.R — 3D CTCRW filter matching supervisor indexing
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
  
  # Observation matrix (3x6): x, y, depth
  Z <- matrix(0, 3, 6)
  Z[1,1] <- 1   # x
  Z[2,3] <- 1   # y
  Z[3,5] <- 1   # depth
  
  # Storage
  aest <- matrix(0, N, 6)
  Pest <- array(0, c(6, 6, N))
  
  v <- matrix(0, N, 3)
  F <- array(0, c(3, 3, N))
  K <- array(0, c(6, 3, N))
  L <- array(0, c(6, 6, N))
  T_array <- array(0, c(6, 6, N))
  
  # NEW: store observation masks
  # This matrix records which observation dimensions were available at every time point.
  obs_mask_store <- matrix(FALSE, N, 3)
  
  # Initial state
  aest[1,] <- a
  Pest[,,1] <- P
  
  ll <- 0
  
  # Loop over time: update at i, predict to i+1
  for (i in 1:(N-1)) {
    
    # Transition matrix
    Tmat <- makeT(
      beta_horiz = beta_horiz_vec[i],
      beta_vert  = beta_vert_vec[i],
      delta      = delta[i+1]
    )
    
    T_array[,,i] <- Tmat
    
    # Process noise matrix
    Qmat <- makeQ(
      beta_horiz = beta_horiz_vec[i],
      beta_vert  = beta_vert_vec[i],
      var_horiz  = var_horiz_vec[i],
      var_vert   = var_vert_vec[i],
      delta      = delta[i+1]
    )
    
    # Observation mask
    obs_mask <- !is.na(y[i,])
    
    # NEW: Store mask for smoother
    obs_mask_store[i,] <- obs_mask
    
    # No observations
    if (!any(obs_mask)) {
      
      aest[i+1,] <- Tmat %*% aest[i,]
      
      Pest[,,i+1] <- Tmat %*% Pest[,,i] %*% t(Tmat) + Qmat
      
      L[,,i] <- Tmat
      
      next
    }
    
    # Observation submatrices
    Z_i <- Z[obs_mask, , drop = FALSE]
    
    H_i <- diag(Hmat[i, obs_mask],sum(obs_mask))
    
    # Innovation
    v_i <- as.numeric(y[i, obs_mask] - Z_i %*% aest[i,])
    
    # Innovation covariance
    F_i <- Z_i %*% Pest[,,i] %*% t(Z_i) + H_i
    
    invF_i <- solve(F_i)
    
    # Kalman gain
    K_i <- Tmat %*% Pest[,,i] %*% t(Z_i) %*% invF_i
    
    # Store innovation and covariance
    v[i, obs_mask] <- v_i
    
    F[obs_mask, obs_mask, i] <- F_i
    
    # Full gain matrix
    K_full <- matrix(0, 6, 3)
    K_full[, obs_mask] <- K_i
    
    K[,,i] <- K_full
    
    # L matrix
    L[,,i] <- Tmat - K_full %*% Z
    
    # Filtered mean
    aest[i+1,] <- Tmat %*% aest[i,] + K_full %*% v[i,]
    
    # Filtered covariance
    Pest[,,i+1] <- Tmat %*% Pest[,,i] %*% t(L[,,i]) + Qmat
    
    # Log-likelihood
    ll <- ll - 0.5 * (log(det(F_i)) + t(v_i) %*% invF_i %*% v_i)
  }
  
  # Return obs_mask in the output
  list(
    ll       = ll,
    aest     = aest,
    Pest     = Pest,
    v        = v,
    F        = F,
    K        = K,
    L        = L,
    T_array  = T_array,
    obs_mask = obs_mask_store
  )
}




















