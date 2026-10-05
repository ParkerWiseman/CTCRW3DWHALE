


############################################################
# CTCRW_smoother.R — RTS smoother matching supervisor indexing
############################################################

CTCRW_smoother <- function(filter_out,
                           beta_horiz_vec,
                           beta_vert_vec,
                           var_horiz_vec,
                           var_vert_vec,
                           delta) {
  
  aest    <- filter_out$aest
  Pest    <- filter_out$Pest
  T_array <- filter_out$T_array
  
  N <- nrow(aest)
  
  a_s <- aest
  P_s <- Pest
  
  for (i in (N-1):1) {
    
    Tmat <- T_array[,,i]
    
    # Recompute predicted covariance (supervisor indexing)
    Qmat <- makeQ(
      beta_horiz        = beta_horiz_vec[i],
      beta_vert        = beta_vert_vec[i],
      var_horiz = var_horiz_vec[i],
      var_vert  = var_vert_vec[i],
      delta     = delta[i+1]
    )
    
    P_pred_next <- Tmat %*% Pest[,,i] %*% t(Tmat) + Qmat
    a_pred_next <- Tmat %*% aest[i,]
    
    # RTS gain
    J <- tryCatch(
      Pest[,,i] %*% t(Tmat) %*% solve(P_pred_next),
      error = function(e) NULL
    )
    
    if (is.null(J)) next
    
    a_s[i,] <- aest[i,] + J %*% (a_s[i+1,] - a_pred_next)
    P_s[,,i] <- Pest[,,i] + J %*% (P_s[,,i+1] - P_pred_next) %*% t(J)
  }
  
  list(a_s = a_s, P_s = P_s)
}








  
  
  
  
  
  
  
  








