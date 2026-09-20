


############################################################
# neg_loglikelihood.R
#
# Negative log-likelihood for the 3D CTCRW model
############################################################


build_Hmat_LinearError <- function(
    data_aug,
    var0_xy,
    var1_xy,
    sd_depth
) {
  
  hd <- data_aug$hdop
  N <- nrow(data_aug)
  
  Hmat <- matrix(
    NA_real_,
    N,
    3
  )
  
  for (i in seq_len(N)) {
    
    h_i <- hd[i]
    
    if (is.na(h_i)) {
      h_i <- 10
    }
    
    var_xy_i <-
      var0_xy +
      var1_xy * h_i
    
    var_xy_i <-
      max(var_xy_i, 1e-3)
    
    Hmat[i,1] <- var_xy_i
    Hmat[i,2] <- var_xy_i
    Hmat[i,3] <- sd_depth^2
  }
  
  Hmat
}


############################################################
# NEGATIVE LOG-LIKELIHOOD
############################################################

neg_loglikelihood <- function(
    params,
    data_aug,
    error_model = c("linearerror")
) {
  
  error_model <- match.arg(error_model)
  
  
  ##########################################################
  # Transform parameters to positive values
  ##########################################################
  
  beta1 <- exp(params["beta1"])
  beta2 <- exp(params["beta2"])
  
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  
  ##########################################################
  # Observations
  ##########################################################
  
  y <- as.matrix(
    data_aug[, c("x", "y", "depth")]
  )
  
  N <- nrow(data_aug)
  
  
  ##########################################################
  # Time in seconds
  ##########################################################
  
  time_sec <-
    as.numeric(
      difftime(
        data_aug$time,
        min(data_aug$time),
        units = "secs"
      )
    )
  
  
  ##########################################################
  # Correctly aligned time intervals
  #
  # First observation = initial state
  # Therefore delta[1] = 0 (first delta is 0, which is correct)
  #
  # No artificial movement at the first observation
  # No forced minimum of 1 second
  # This change alone makes the likelihood smoother and more stable
  ##########################################################
  
  delta <- c(0,diff(time_sec))
  
  delta[!is.finite(delta)] <- 0
  
  delta <- pmax(delta, 0)
  
  
  ##########################################################
  # Process-noise variances
  ##########################################################
  
  s_horiz <- sigma1^2
  s_vert <- sigma2^2
  
  ##########################################################
  # Parameter vectors
  ##########################################################
  
  beta1_vec <- rep(beta1,N)
  
  beta2_vec <- rep(beta2,N)
  
  ##########################################################
  # Measurement error model
  ##########################################################
  
  if (error_model == "linearerror") {
    
    Hmat <- build_Hmat_LinearError(
      data_aug = data_aug,
      var0_xy = 100,
      var1_xy = 2,
      sd_depth = 25
    )
  }
  
  ##########################################################
  # Initial state
  ##########################################################
  
  get_first_non_missing <- function(x) {
    
    idx <- which(!is.na(x))[1]
    
    if (length(idx) == 0) {
      return(0)
    }
    
    x[idx]
  }
  
  
  a <- c(
    get_first_non_missing(y[,1]),
    0,
    
    get_first_non_missing(y[,2]),
    0,
    
    get_first_non_missing(y[,3]),
    0
  )
  
  ##########################################################
  # Initial covariance
  ##########################################################
  
  P <- diag(6) * 1e6
  
  ##########################################################
  # Kalman filter
  ##########################################################
  
  filt <- CTCRW_filter1(
    y = y,
    Hmat = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz = s_horiz,
    s_vert = s_vert,
    delta = delta,
    a = a,
    P = P
  )
  
  ##########################################################
  # Negative log-likelihood
  ##########################################################
  
  if (!is.finite(filt$ll)) {
    return(1e100)
  }
  
  -filt$ll
}













