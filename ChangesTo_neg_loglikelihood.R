





neg_loglikelihood <- function(params, data_aug, error_model = c("noerror","constanterror","linearerror")) {
  
  error_model <- match.arg(error_model)
  
  # unpack parameters
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  # data
  y <- as.matrix(data_aug[, c("x","y","depth")])
  N <- nrow(data_aug)
  
  # ---------------------------------------------------------
  # STABILIZED DELTA (CRITICAL FIX)
  # ---------------------------------------------------------
  delta_raw   <- diff(data_aug$Time)
  delta_fixed <- pmax(delta_raw, 1e-5)   # enforce minimum positive time step
  delta       <- c(delta_fixed[1], delta_fixed)
  
  # process noise (2-sigma model)
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, N)
  beta2_vec <- rep(beta2, N)
  
  # ---------------------------------------------------------
  # BUILD H MATRIX BASED ON error_model
  # ---------------------------------------------------------
  
  if (error_model == "noerror") {
    
    eps <- 1e-6
    Hmat <- matrix(eps, N, 3)
    
  } else if (error_model == "constanterror") {
    
    sd_xy    <- 5
    sd_depth <- 10
    
    Hmat <- matrix(NA_real_, N, 3)
    Hmat[,1] <- sd_xy^2
    Hmat[,2] <- sd_xy^2
    Hmat[,3] <- sd_depth^2
    
  } else if (error_model == "linearerror") {
    
    var0_xy  <- 0
    var1_xy  <- 0.1
    sd_depth <- 10
    
    hd <- data_aug$hdop
    
    Hmat <- matrix(NA_real_, N, 3)
    
    for (i in 1:N) {
      if (!is.na(hd[i])) {
        var_xy_i <- var0_xy + var1_xy * hd[i]
        Hmat[i,1] <- var_xy_i
        Hmat[i,2] <- var_xy_i
      } else {
        Hmat[i,1] <- 1e6
        Hmat[i,2] <- 1e6
      }
      Hmat[i,3] <- sd_depth^2
    }
  }
  
  # ---------------------------------------------------------
  # INITIAL STATE
  # ---------------------------------------------------------
  
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  
  P <- diag(6) * 1e2
  
  # ---------------------------------------------------------
  # FILTER CALL
  # ---------------------------------------------------------
  
  filt <- CTCRW_filter1(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz   = s_horiz,
    s_vert    = s_vert,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  # ---------------------------------------------------------
  # SAFETY CHECK — PREVENT NaN/Inf FROM CRASHING optim()
  # ---------------------------------------------------------
  
  if (!is.finite(filt$ll)) {
    return(1e10)
  }
  
  -filt$ll
}





###################################################


neg_loglikelihood <- function(params, data_aug, 
                              error_model = c("noerror","constanterror","linearerror")) {
  
  error_model <- match.arg(error_model)
  
  # unpack parameters
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  # data
  y <- as.matrix(data_aug[, c("x","y","depth")])
  N <- nrow(data_aug)
  
  # STABILIZED DELTA
  delta_raw   <- diff(data_aug$Time)
  delta_fixed <- pmax(delta_raw, 1e-5)
  delta       <- c(delta_fixed[1], delta_fixed)
  
  # process noise
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, N)
  beta2_vec <- rep(beta2, N)
  
  # BUILD H MATRIX
  if (error_model == "noerror") {
    
    eps <- 1e-6
    Hmat <- matrix(eps, N, 3)
    
  } else if (error_model == "constanterror") {
    
    sd_xy    <- 5
    sd_depth <- 10
    
    Hmat <- matrix(NA_real_, N, 3)
    Hmat[,1] <- sd_xy^2
    Hmat[,2] <- sd_xy^2
    Hmat[,3] <- sd_depth^2
    
  } else if (error_model == "linearerror") {
    
    var0_xy  <- 0
    var1_xy  <- 0.1
    sd_depth <- 10
    
    hd <- data_aug$hdop
    
    Hmat <- matrix(NA_real_, N, 3)
    
    for (i in 1:N) {
      if (!is.na(hd[i])) {
        var_xy_i <- var0_xy + var1_xy * hd[i]
        Hmat[i,1] <- var_xy_i
        Hmat[i,2] <- var_xy_i
      } else {
        Hmat[i,1] <- 1e6
        Hmat[i,2] <- 1e6
      }
      Hmat[i,3] <- sd_depth^2
    }
  }
  
  # INITIAL STATE
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  
  P <- diag(6) * 1e2
  
  # FILTER
  filt <- CTCRW_filter1(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz   = s_horiz,
    s_vert    = s_vert,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  # BASE NLL
  nll <- -filt$ll
  
  # SOFT PENALTY (no hard bounds)
  # penalize absurdly large or tiny values, but smoothly
  pen <- 0
  
  # center on reasonable scales, penalize log-distance
  pen <- pen + 0.01 * (max(0, log(beta1) - log(10))^2 + max(0, log(0.01) - log(beta1))^2)
  pen <- pen + 0.01 * (max(0, log(beta2) - log(10))^2 + max(0, log(0.01) - log(beta2))^2)
  pen <- pen + 0.01 * (max(0, log(sigma1) - log(200))^2 + max(0, log(0.1) - log(sigma1))^2)
  pen <- pen + 0.01 * (max(0, log(sigma2) - log(200))^2 + max(0, log(0.1) - log(sigma2))^2)
  
  # guard against non-finite LL
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  nll + pen
}



#################################





neg_loglikelihood_whale_beta_only <- function(params, data_aug,
                                              sigma1_fixed = 30,
                                              sigma2_fixed = 10,
                                              error_model = c("linearerror")) {
  
  error_model <- match.arg(error_model)
  
  # unpack parameters: only beta1, beta2 are free
  beta1 <- exp(params["beta1"])
  beta2 <- exp(params["beta2"])
  
  # fixed process noise from prior knowledge / simulations
  sigma1 <- sigma1_fixed
  sigma2 <- sigma2_fixed
  
  # data
  y <- as.matrix(data_aug[, c("x","y","depth")])
  N <- nrow(data_aug)
  
  # stabilized delta
  delta_raw   <- diff(data_aug$Time)
  delta_fixed <- pmax(delta_raw, 1e-5)
  delta       <- c(delta_fixed[1], delta_fixed)
  
  # process noise (2-sigma model)
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, N)
  beta2_vec <- rep(beta2, N)
  
  # H matrix for linear HDOP
  var0_xy  <- 0
  var1_xy  <- 0.1
  sd_depth <- 10
  
  hd <- data_aug$hdop
  
  Hmat <- matrix(NA_real_, N, 3)
  for (i in 1:N) {
    if (!is.na(hd[i])) {
      var_xy_i <- var0_xy + var1_xy * hd[i]
      Hmat[i,1] <- var_xy_i
      Hmat[i,2] <- var_xy_i
    } else {
      Hmat[i,1] <- 1e6
      Hmat[i,2] <- 1e6
    }
    Hmat[i,3] <- sd_depth^2
  }
  
  # initial state
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  P <- diag(6) * 1e2
  
  # filter
  filt <- CTCRW_filter1(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz   = s_horiz,
    s_vert    = s_vert,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  nll <- -filt$ll
  
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  nll
}























###################################################







###########################
# DELETE BELOW CODE:


neg_loglikelihood <- function(params, data_aug,
                              error_model = c("noerror","constanterror","linearerror")) {
  
  error_model <- match.arg(error_model)
  
  # unpack parameters
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  # data
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  # ---- SAFE DELTA ----
  delta_raw <- c(NA, diff(data_aug$Time))
  # replace first NA with second value
  delta_raw[1] <- delta_raw[2]
  # enforce strictly positive, finite time steps
  delta <- pmax(delta_raw, 1e-6)
  delta[!is.finite(delta)] <- 1e-6
  
  # process noise (2-sigma model)
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### ----------------------------------------------------
  #### BUILD H MATRIX BASED ON error_model (ROBUST)
  #### ----------------------------------------------------
  
  if (error_model == "noerror") {
    
    eps <- 1e-6
    Hmat <- matrix(eps, nrow(data_aug), 3)
    
  } else if (error_model == "constanterror") {
    
    sd_xy    <- 5
    sd_depth <- 10
    
    Hmat <- matrix(sd_xy^2, nrow(data_aug), 3)
    Hmat[,3] <- sd_depth^2
    
  } else if (error_model == "linearerror") {
    
    var0_xy  <- 0
    var1_xy  <- 0.1
    sd_depth <- 10
    
    hd <- data_aug$hdop
    N  <- nrow(data_aug)
    
    Hmat <- matrix(NA_real_, N, 3)
    
    for (i in 1:N) {
      # clean HDOP: replace NA with large value, enforce minimum
      h_i <- hd[i]
      if (is.na(h_i)) {
        h_i <- 50    # very uncertain GPS
      }
      h_i <- max(h_i, 1)  # avoid zero/negative HDOP
      
      var_xy_i <- var0_xy + var1_xy * h_i
      var_xy_i <- max(var_xy_i, 1e-6)  # avoid zero/negative variance
      
      Hmat[i,1] <- var_xy_i
      Hmat[i,2] <- var_xy_i
      Hmat[i,3] <- sd_depth^2
    }
  }
  
  #### ----------------------------------------------------
  #### INITIAL STATE
  #### ----------------------------------------------------
  
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  
  P <- diag(6) * 1e2
  
  #### ----------------------------------------------------
  #### FILTER CALL
  #### ----------------------------------------------------
  
  filt <- CTCRW_filter1(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz   = s_horiz,
    s_vert    = s_vert,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  #### ----------------------------------------------------
  #### RETURN NEGATIVE LOG-LIKELIHOOD
  #### ----------------------------------------------------
  
  -filt$ll
}














############

# AUGUST 20 CODE


neg_loglikelihood <- function(params, data_aug,
                              error_model = c("noerror","constanterror","linearerror")) {
  
  error_model <- match.arg(error_model)
  
  # unpack parameters
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  # data
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  # ---- SAFE DELTA ----
  delta_raw <- c(NA, diff(data_aug$Time))
  delta_raw[1] <- delta_raw[2]
  delta <- pmax(delta_raw, 1e-6)
  delta[!is.finite(delta)] <- 1e-6
  
  # process noise (2-sigma model)
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### ----------------------------------------------------
  #### BUILD H MATRIX BASED ON error_model (ROBUST)
  #### ----------------------------------------------------
  
  if (error_model == "noerror") {
    
    eps <- 1e-6
    Hmat <- matrix(eps, nrow(data_aug), 3)
    
  } else if (error_model == "constanterror") {
    
    sd_xy    <- 5
    sd_depth <- 10
    
    Hmat <- matrix(sd_xy^2, nrow(data_aug), 3)
    Hmat[,3] <- sd_depth^2
    
  } else if (error_model == "linearerror") {
    
    var0_xy  <- 0
    var1_xy  <- 0.1
    sd_depth <- 10
    
    hd <- data_aug$hdop
    N  <- nrow(data_aug)
    
    Hmat <- matrix(NA_real_, N, 3)
    
    for (i in 1:N) {
      # clean HDOP: replace NA with large value, enforce minimum
      h_i <- hd[i]
      if (is.na(h_i)) h_i <- 50
      h_i <- max(h_i, 1)
      
      var_xy_i <- var0_xy + var1_xy * h_i
      var_xy_i <- max(var_xy_i, 1e-6)
      
      Hmat[i,1] <- var_xy_i
      Hmat[i,2] <- var_xy_i
      Hmat[i,3] <- sd_depth^2
    }
  }
  
  #### ----------------------------------------------------
  #### INITIAL STATE
  #### ----------------------------------------------------
  
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  
  P <- diag(6) * 1e2
  
  #### ----------------------------------------------------
  #### FILTER CALL
  #### ----------------------------------------------------
  
  filt <- CTCRW_filter1(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz   = s_horiz,
    s_vert    = s_vert,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  #### ----------------------------------------------------
  #### RETURN NEGATIVE LOG-LIKELIHOOD
  #### ----------------------------------------------------
  
  -filt$ll
}

















