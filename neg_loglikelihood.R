


############################################################
# neg_loglikelihood.R — CTCRW likelihood with delta in seconds
############################################################

build_Hmat_LinearError <- function(data_aug, var0_xy, var1_xy, sd_depth) {
  
  hd <- data_aug$hdop
  N  <- nrow(data_aug)
  
  Hmat <- matrix(NA_real_, N, 3)
  
  for (i in 1:N) {
    
    h_i <- hd[i]
    if (is.na(h_i)) h_i <- 10
    
    var_xy_i <- var0_xy + var1_xy * h_i
    var_xy_i <- max(var_xy_i, 1e-3)
    
    Hmat[i,1] <- var_xy_i
    Hmat[i,2] <- var_xy_i
    Hmat[i,3] <- sd_depth^2
  }
  
  Hmat
}

neg_loglikelihood <- function(params, data_aug,
                              error_model = c("linearerror")) {
  
  error_model <- match.arg(error_model)
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  #### ----------------------------------------------------
  #### DELTA IN SECONDS (minimum = 1 second)
  #### ----------------------------------------------------
  
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  
  delta_raw <- c(NA, diff(time_sec))
  delta_raw[1] <- delta_raw[2]
  
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  #### ----------------------------------------------------
  #### PROCESS NOISE
  #### ----------------------------------------------------
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### ----------------------------------------------------
  #### H MATRIX (your chosen values)
  #### ----------------------------------------------------
  
  var0_xy  <- 100
  var1_xy  <- 2
  sd_depth <- 25
  
  Hmat <- build_Hmat_LinearError(data_aug, var0_xy, var1_xy, sd_depth)
  
  #### ----------------------------------------------------
  #### INITIAL STATE (first non-missing)
  #### ----------------------------------------------------
  
  get_first_non_missing <- function(col) {
    idx <- which(!is.na(col))[1]
    if (is.na(idx)) 0 else col[idx]
  }
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  #### ----------------------------------------------------
  #### INITIAL COVARIANCE
  #### ----------------------------------------------------
  
  P <- diag(6) * 1e6
  
  #### ----------------------------------------------------
  #### FILTER
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
  
  -filt$ll
}









