


############################################################
# neg_loglikelihood.R
# Unified, robust CTCRW negative log-likelihood
# Supports: "noerror", "constanterror", "linearerror"
############################################################

#### ----------------------------------------------------
#### H-MATRIX BUILDERS
#### ----------------------------------------------------

# 1. No error (almost perfect observations)
build_Hmat_NoError <- function(data_aug) {
  N   <- nrow(data_aug)
  eps <- 1e-6  # small variance for numerical stability
  
  Hmat <- matrix(0, N, 3)
  Hmat[,1] <- eps
  Hmat[,2] <- eps
  Hmat[,3] <- eps
  
  Hmat
}

# 2. Constant error (same variance for all points)
build_Hmat_ConstantError <- function(data_aug, sd_xy, sd_depth) {
  N <- nrow(data_aug)
  
  Hmat <- matrix(sd_xy^2, N, 3)
  Hmat[,3] <- sd_depth^2
  
  Hmat
}

# 3. Linear HDOP error (variance grows with HDOP)
build_Hmat_LinearError <- function(data_aug, var0_xy, var1_xy, sd_depth) {
  hd <- data_aug$hdop
  N  <- nrow(data_aug)
  
  Hmat <- matrix(NA_real_, N, 3)
  
  for (i in 1:N) {
    # clean HDOP: replace NA with large value, enforce minimum
    h_i <- hd[i]
    if (is.na(h_i)) h_i <- 50        # very uncertain GPS
    h_i <- max(h_i, 1)               # avoid zero/negative HDOP
    
    var_xy_i <- var0_xy + var1_xy * h_i
    var_xy_i <- max(var_xy_i, 1e-6)  # avoid zero/negative variance
    
    Hmat[i,1] <- var_xy_i
    Hmat[i,2] <- var_xy_i
    Hmat[i,3] <- sd_depth^2
  }
  
  Hmat
}

#### ----------------------------------------------------
#### ROBUST UNIFIED NEGATIVE LOG-LIKELIHOOD (AUG 20)
#### ----------------------------------------------------

neg_loglikelihood <- function(params, data_aug,
                              error_model = c("noerror","constanterror","linearerror")) {
  
  error_model <- match.arg(error_model)
  
  # unpack parameters (log-scale → positive)
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])   # horizontal process SD
  sigma2 <- exp(params["sigma2"])   # vertical process SD
  
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
    
    Hmat <- build_Hmat_NoError(data_aug)
    
  } else if (error_model == "constanterror") {
    
    sd_xy    <- 5    # e.g., meters
    sd_depth <- 10   # e.g., meters
    
    Hmat <- build_Hmat_ConstantError(data_aug, sd_xy, sd_depth)
    
  } else if (error_model == "linearerror") {
    
    var0_xy  <- 0    # baseline variance
    var1_xy  <- 0.1  # variance per HDOP unit
    sd_depth <- 10   # depth SD
    
    Hmat <- build_Hmat_LinearError(data_aug, var0_xy, var1_xy, sd_depth)
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















