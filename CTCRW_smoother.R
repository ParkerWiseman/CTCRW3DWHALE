


############################################################
# CTCRW_smoother.R
# Rauch-Tung-Striebel smoother for the 3D CTCRW model
# State: [x, vx, y, vy, depth, vdepth]
############################################################

CTCRW_smoother1 <- function(filter_out, beta1_vec, beta2_vec, s_horiz, s_vert, delta) {
  N <- nrow(filter_out$a_f)
  
  ##########################################################
  # Start with filtered estimates
  ##########################################################
  
  a_s <- filter_out$a_f
  P_s <- filter_out$P_f
  
  ##########################################################
  # Backward smoothing
  ##########################################################
  
  if (N >= 2) {
    
    for (i in (N - 1):1) {
      
      ######################################################
      # Transition from i to i+1
      ######################################################
      
      Tmat <- makeT(b1 = beta1_vec[i + 1],b2 = beta2_vec[i + 1],delta = delta[i + 1])
      
      ######################################################
      # Filtered state at time i
      ######################################################
      
      a_f_i <- filter_out$a_f[i, ]
      
      ######################################################
      # Predicted state at time i+1
      ######################################################
      
      a_p_ip1 <- filter_out$a_p[i + 1, ]
      
      ######################################################
      # Covariances
      ######################################################
      
      P_f_i <- filter_out$P_f[[i]]
      P_p_ip1 <- filter_out$P_p[[i + 1]]
      
      ######################################################
      # Numerical stabilization
      # Symmetrizing covariances:
      # prevents asymmetry drift from floating‑point accumulation
      ######################################################
      
      P_p_ip1 <- (P_p_ip1 + t(P_p_ip1)) / 2
      P_f_i <- (P_f_i + t(P_f_i)) / 2
      
      ######################################################
      # Smoothing gain
      # J_i = P_f_i T' P_p_(i+1)^(-1)
      # Adding jitter before inversion:
      # prevents singularities and near‑singular matrices from blowing up the smoother
      ######################################################
      
      P_p_jitter <- P_p_ip1 + diag(1e-8, 6)
      J <- tryCatch(P_f_i %*% t(Tmat) %*% solve(P_p_jitter), error = function(e) NULL)
      
      ######################################################
      # Checking for numerical failure:
      # If numerical failure occurs, retain filtered state
      ######################################################
      
      if (is.null(J)) {
        next
      }
      
      ######################################################
      # Smoothed state
      ######################################################
      
      a_s[i, ] <- as.numeric(a_f_i + J %*% (a_s[i + 1, ] - a_p_ip1))
      
      ######################################################
      # Smoothed covariance
      # Covariance symmetrization after smoothing
      ######################################################
      
      P_s[[i]] <- P_f_i + J %*% (P_s[[i + 1]] - P_p_ip1) %*% t(J)
      P_s[[i]] <- (P_s[[i]] + t(P_s[[i]])) / 2
    }
  }
  
  ##########################################################
  # Return
  ##########################################################
  
  list(a_s = a_s,P_s = P_s)
}

############################################################
# Backwards-compatible version
############################################################

CTCRW_smoother <- function(filter_out,beta1_vec,beta2_vec,s1,s2,s3,delta) {
  
  CTCRW_smoother1(filter_out = filter_out,beta1_vec = beta1_vec,beta2_vec = beta2_vec,s_horiz = s1,s_vert = s3,delta = delta)
}

  
  
  
  
  
  
  
  








