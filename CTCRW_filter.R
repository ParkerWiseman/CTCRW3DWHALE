


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



























