


############################################################
# CTCRW_smoother.R
#
# Durbin & Koopman disturbance smoother
# 3D CTCRW version with partial-observation handling
#
# Updated so that the final observation at time N
# contributes to the smoothing recursion.
############################################################

CTCRW_smoother <- function(filter_out) {
  
  aest <- filter_out$aest
  Pest <- filter_out$Pest
  v    <- filter_out$v
  F    <- filter_out$F
  L    <- filter_out$L
  
  obs_mask_store <- filter_out$obs_mask
  
  Nobs <- nrow(aest)
  
  ##########################################################
  # Observation matrix
  ##########################################################
  
  Z <- matrix(0, 3, 6)
  
  Z[1,1] <- 1   # x
  Z[2,3] <- 1   # y
  Z[3,5] <- 1   # depth
  
  ##########################################################
  # Disturbance smoother variables
  ##########################################################
  
  r <- matrix(0, 6, 1)
  
  Nmat <- matrix(0, 6, 6)
  
  ##########################################################
  # Incorporate final observation at time N
  ##########################################################
  
  obs_mask_N <- obs_mask_store[Nobs, ]
  
  if (any(obs_mask_N)) {
    
    Z_use <- Z[obs_mask_N, , drop = FALSE]
    
    v_use <- matrix(v[Nobs, obs_mask_N], ncol = 1)
    
    F_use <- F[obs_mask_N,obs_mask_N,Nobs,drop = FALSE]
    
    Finv <- solve(F_use)
    
    r <- t(Z_use) %*% Finv %*% v_use
    
    Nmat <- t(Z_use) %*% Finv %*% Z_use
  }
  
  ##########################################################
  # Output
  ##########################################################
  
  a_s <- aest
  P_s <- Pest
  
  ##########################################################
  # Backward recursion
  ##########################################################
  
  for (i in Nobs:2) {
    
    obs_mask <- obs_mask_store[i - 1, ]
    
    ########################################################
    # No observations available
    ########################################################
    
    if (!any(obs_mask)) {
      
      r <- t(L[,,i-1]) %*% r
      
      Nmat <- t(L[,,i-1]) %*% Nmat %*% L[,,i-1]
      
    } else {
      
      ######################################################
      # Observed dimensions only
      ######################################################
      
      v_use <- matrix(v[i-1, obs_mask],ncol = 1)
      
      Z_use <- Z[obs_mask,,drop = FALSE]
      
      F_full <- F[,,i-1]
      
      F_use <- F_full[obs_mask,obs_mask,drop = FALSE]
      
      Finv <- solve(F_use)
      
      ######################################################
      # Disturbance smoother update
      ######################################################
      
      r <- t(Z_use) %*% (Finv %*% v_use) + t(L[,,i-1]) %*% r
      
      Nmat <- t(Z_use) %*% Finv %*% Z_use + t(L[,,i-1]) %*% Nmat %*% L[,,i-1]
    }
    
    ########################################################
    # Smoothed moments
    ########################################################
    
    a_s[i-1, ] <- aest[i-1, ] + as.vector(Pest[,,i-1] %*% r)
    
    P_s[,,i-1] <- Pest[,,i-1] - Pest[,,i-1] %*% Nmat %*% Pest[,,i-1]
  }
  
  list(
    a_s = a_s,
    P_s = P_s
  )
}






