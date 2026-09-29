


############################################################
# CTCRW_smoother.R
############################################################

CTCRW_smoother1 <- function(
    filter_out,
    beta1_vec,
    beta2_vec,
    s_horiz,
    s_vert,
    delta
) {
  
  N <- nrow(filter_out$a_f)
  
  a_s <- filter_out$a_f
  P_s <- filter_out$P_f
  
  # Use stored T matrices
  T_array <- filter_out$T_array
  
  if (N >= 2) {
    
    for (i in (N - 1):1) {
      
      # Use stored Tmat
      Tmat <- T_array[[i + 1]]
      
      a_f_i   <- filter_out$a_f[i, ]
      a_p_ip1 <- filter_out$a_p[i + 1, ]
      
      P_f_i   <- filter_out$P_f[[i]]
      P_p_ip1 <- filter_out$P_p[[i + 1]]
      
      # RTS smoothing gain
      J <- tryCatch(
        P_f_i %*% t(Tmat) %*% solve(P_p_ip1),
        error = function(e) NULL
      )
      
      if (is.null(J)) next
      
      # Smoothed mean
      a_s[i, ] <- as.numeric(
        a_f_i + J %*% (a_s[i + 1, ] - a_p_ip1)
      )
      
      # Smoothed covariance
      P_s[[i]] <- P_f_i + J %*% (P_s[[i + 1]] - P_p_ip1) %*% t(J)
    }
  }
  
  list(a_s = a_s, P_s = P_s)
}

CTCRW_smoother <- function(
    filter_out,
    beta1_vec,
    beta2_vec,
    s1,
    s2,
    s3,
    delta
) {
  CTCRW_smoother1(
    filter_out = filter_out,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz = s1,
    s_vert = s3,
    delta = delta
  )
}











  
  
  
  
  
  
  
  








