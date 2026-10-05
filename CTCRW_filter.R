


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
  aest <- matrix(0, N, 6)              # filtered means
  Pest <- array(0, c(6, 6, N))         # filtered covariances
  
  v <- matrix(0, N, 3)                 # innovations (full 3D)
  F <- array(0, c(3, 3, N))            # innovation covariances
  K <- array(0, c(6, 3, N))            # Kalman gains
  L <- array(0, c(6, 6, N))            # L matrices
  T_array <- array(0, c(6, 6, N))      # transition matrices
  
  # Initial state
  aest[1,] <- a
  Pest[,,1] <- P
  
  ll <- 0
  
  # Loop over time: update at time i, predict to time i+1
  for (i in 1:(N-1)) {
    
    # Transition matrix for step i -> i+1
    Tmat <- makeT(beta_horiz = beta_horiz_vec[i], beta_vert = beta_vert_vec[i], delta = delta[i+1])
    
    T_array[,,i] <- Tmat
    
    # Process noise for step i -> i+1
    Qmat <- makeQ(beta_horiz = beta_horiz_vec[i], beta_vert = beta_vert_vec[i],
      var_horiz = var_horiz_vec[i], var_vert = var_vert_vec[i], delta = delta[i+1]
    )
    
    # Missing observation at time i?
    obs_mask <- !is.na(y[i,])
    
    #if (!any(obs_mask)) {
    #  # No observation: pure prediction
    #  aest[i+1,]   <- Tmat %*% aest[i,]
    #  Pest[,,i+1]  <- Tmat %*% Pest[,,i] %*% t(Tmat) + Qmat
    #  L[,,i]       <- Tmat
    #  next
    #}
    
    # Subset Z and H to observed components at time i
    Z_i <- Z[obs_mask,, drop = FALSE]
    H_i <- diag(Hmat[i, obs_mask], sum(obs_mask))
    
    # Innovation only for observed components (update at time i)
    v_i <- as.numeric(y[i, obs_mask] - Z_i %*% aest[i,])
    
    # Innovation covariance
    F_i <- Z_i %*% Pest[,,i] %*% t(Z_i) + H_i
    invF_i <- solve(F_i)
    
    # Kalman gain for observed components
    K_i <- Tmat %*% Pest[,,i] %*% t(Z_i) %*% invF_i  # 6 x k
    
    # Embed v_i and F_i back into full 3D space at time i
    v[i, obs_mask] <- v_i
    F[obs_mask, obs_mask, i] <- F_i
    
    # Full 6x3 K: zeros then fill observed columns
    K_full <- matrix(0, 6, 3)
    K_full[, obs_mask] <- K_i
    K[,,i] <- K_full
    
    # L matrix for smoother
    L[,,i] <- Tmat - K_full %*% Z
    
    # Filtered mean at time i+1
    aest[i+1,] <- Tmat %*% aest[i,] + K_full %*% v[i,]
    
    # Filtered covariance at time i+1 (DK / crawl form)
    Pest[,,i+1] <- Tmat %*% Pest[,,i] %*% t(L[,,i]) + Qmat
    
    # Log-likelihood contribution (only observed components)
    ll <- ll - 0.5 * (log(det(F_i)) + t(v_i) %*% invF_i %*% v_i)
  }
  
  list(
    ll      = ll,
    aest    = aest,
    Pest    = Pest,
    v       = v,
    F       = F,
    K       = K,
    L       = L,
    T_array = T_array
  )
}
























