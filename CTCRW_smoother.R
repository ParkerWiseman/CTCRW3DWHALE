CTCRW_smoother <- function(filter_out, beta1_vec, beta2_vec,
                           s1, s2, s3, delta) {
  
  N <- nrow(filter_out$a_f)
  a_s <- filter_out$a_f
  P_s <- filter_out$P_f
  
  for (i in (N-1):1) {
    
    Tmat <- makeT_R(beta1_vec[i+1], beta2_vec[i+1], delta[i+1])
    
    a_f_i   <- filter_out$a_f[i,]
    a_p_ip1 <- filter_out$a_p[i+1,]
    P_f_i   <- filter_out$P_f[[i]]
    P_p_ip1 <- filter_out$P_p[[i+1]]
    
    J <- tryCatch(P_f_i %*% t(Tmat) %*% solve(P_p_ip1), error=function(e) NULL)
    if (is.null(J)) next
    
    a_s[i,]  <- a_f_i + J %*% (a_s[i+1,] - a_p_ip1)
    P_s[[i]] <- P_f_i + J %*% (P_s[[i+1]] - P_p_ip1) %*% t(J)
  }
  
  list(a_s = a_s, P_s = P_s)
}