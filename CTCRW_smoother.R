


############################################################
# CTCRW_smoother.R — Durbin & Koopman disturbance smoother
# Matching supervisor indexing, extended to 6‑dimensional 3D CTCRW
############################################################

CTCRW_smoother <- function(filter_out) {
  
  aest <- filter_out$aest      # filtered means (N × 6)
  Pest <- filter_out$Pest      # filtered covariances (6×6×N)
  v    <- filter_out$v         # innovations (N × 3)
  F    <- filter_out$F         # innovation covariances (3×3×N)
  L    <- filter_out$L         # L matrices (6×6×N)
  
  Nobs <- nrow(aest)
  
  # Observation matrix (same as filter)
  Z <- matrix(0, 3, 6)
  Z[1,1] <- 1   # x
  Z[2,3] <- 1   # y
  Z[3,5] <- 1   # depth
  
  # Disturbance smoother variables
  r <- rep(0, 6)
  Nmat <- matrix(0, 6, 6)
  
  # Output smoothed means/covariances
  a_s <- aest
  P_s <- Pest
  
  # Backward loop (Durbin & Koopman)
  for (i in Nobs:2) {
    
    # Missing observation check
    missing_i <- any(is.na(v[i-1,])) || any(F[,,i-1] == 0)
    
    if (missing_i) {
      
      # No update from observation
      r    <- t(L[,,i-1]) %*% r
      Nmat <- t(L[,,i-1]) %*% Nmat %*% L[,,i-1]
      
    } else {
      
      # Disturbance smoother update
      Finv <- solve(F[,,i-1])
      
      r <- t(Z) %*% (Finv %*% v[i-1,]) + t(L[,,i-1]) %*% r
      Nmat <- t(Z) %*% Finv %*% Z + t(L[,,i-1]) %*% Nmat %*% L[,,i-1]
    }
    
    # Smoothed mean and covariance
    a_s[i-1,] <- aest[i-1,] + Pest[,,i-1] %*% r
    P_s[,,i-1] <- Pest[,,i-1] - Pest[,,i-1] %*% Nmat %*% Pest[,,i-1]
  }
  
  list(a_s = a_s, P_s = P_s)
}





  
  
  
  
  
  
  
  








