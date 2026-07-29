

#### 1. NO ERROR SITUATION ####




build_Hmat_NoError <- function(data_aug) {
  N <- nrow(data_aug)
  
  # Perfect observations → extremely small variance (not zero)
  # Zero variance breaks the Kalman filter (det(F) = 0)
  eps <- 1e-6
  
  Hmat <- matrix(0, N, 3)
  Hmat[,1] <- eps
  Hmat[,2] <- eps
  Hmat[,3] <- eps
  
  Hmat
}




neg_loglikelihood_NoError1 <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  delta <- rep(dt, nrow(data_aug))
  
  s1 <- sigma1^2
  s2 <- sigma1^2
  s3 <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  a <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
  P <- diag(6) * 1e2
  
  filt <- CTCRW_filter(
    y         = y,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s1        = s1,
    s2        = s2,
    s3        = s3,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  -filt$ll
}



















############################################







neg_loglikelihood <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  
  var0_xy  <- 5^2
  var1_xy  <- 1^2
  sd_depth <- 10
  
  Hmat <- build_Hmat(data_aug, var0_xy, var1_xy, sd_depth)
  
  delta <- rep(dt, nrow(data_aug))
  
  s1 <- sigma1^2
  s2 <- sigma1^2
  s3 <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  P <- diag(6) * 1e2
  
  filt <- CTCRW_filter(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s1        = s1,
    s2        = s2,
    s3        = s3,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  -filt$ll
}



###########################################################

#### 2. CONSTANT ERROR SITUATION ####


#build_Hmat_ConstantError <- function(N, sd_xy, sd_depth) {
#  Hmat <- matrix(NA_real_, N, 3)
#  Hmat[,1] <- sd_xy^2
#  Hmat[,2] <- sd_xy^2
#  Hmat[,3] <- sd_depth^2
#  Hmat
#}

build_Hmat_ConstantError <- function(data_aug, var0_xy, var1_xy, sd_depth) {
  N <- nrow(data_aug)
  Hmat <- matrix(NA_real_, N, 3)
  Hmat[,1] <- var0_xy
  Hmat[,2] <- var0_xy
  Hmat[,3] <- sd_depth^2
  Hmat
}




neg_loglikelihood_ConstantError <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  
  var0_xy  <- 5^2
  var1_xy  <- 1^2
  sd_depth <- 10
  
  Hmat <- build_Hmat_ConstantError(data_aug, var0_xy, var1_xy, sd_depth)
  
  delta <- rep(dt, nrow(data_aug))
  
  s1 <- sigma1^2
  s2 <- sigma1^2
  s3 <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  P <- diag(6) * 1e2
  
  filt <- CTCRW_filter(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s1        = s1,
    s2        = s2,
    s3        = s3,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  -filt$ll
}


# (using) s_horiz and s_vert) #

build_Hmat_ConstantError2 <- function(data_aug, sd_xy, sd_depth) {
  N <- nrow(data_aug)
  Hmat <- matrix(NA_real_, N, 3)
  Hmat[,1] <- sd_xy^2
  Hmat[,2] <- sd_xy^2
  Hmat[,3] <- sd_depth^2
  Hmat
}


# (using) s_horiz and s_vert) #

neg_loglikelihood_ConstantError2 <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])   # horizontal
  sigma2 <- exp(params["sigma2"])   # vertical
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  delta <- rep(dt, nrow(data_aug))
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  Hmat <- build_Hmat_ConstantError2(data_aug, sd_xy, sd_depth)
  
  a <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
  P <- diag(6) * 1e2
  
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









############################################################


# LINEAR ERROR SITUATION #


build_Hmat_LinearError <- function(data_aug, var0_xy, var1_xy, sd_depth) {
  y  <- as.matrix(data_aug[, c("x","y","depth")])
  hd <- data_aug$hdop
  
  Hmat <- matrix(NA_real_, nrow(data_aug), 3)
  
  for (i in 1:nrow(data_aug)) {
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
  
  Hmat
}














##################################
##################################
##################################
##################################

# (functions using s_horiz and s_vert)




#(same as build_Hmat_NoError)

build_Hmat_NoError1 <- function(data_aug) {
  N <- nrow(data_aug)
  
  eps <- 1e-6   # small variance for numerical stability
  
  Hmat <- matrix(0, N, 3)
  Hmat[,1] <- eps
  Hmat[,2] <- eps
  Hmat[,3] <- eps
  
  Hmat
}




neg_loglikelihood_NoError1 <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])   # horizontal
  sigma2 <- exp(params["sigma2"])   # vertical
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  delta <- rep(dt, nrow(data_aug))
  
  # NEW: only two process-noise variances
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  a <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
  P <- diag(6) * 1e2
  
  Hmat <- build_Hmat_NoError(data_aug)
  
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









# (Uses s_horiz and s_vert)

build_Hmat_LinearError2 <- function(data_aug, var0_xy, var1_xy, sd_depth) {
  
  hd <- data_aug$hdop
  N  <- nrow(data_aug)
  
  Hmat <- matrix(NA_real_, N, 3)
  
  for (i in 1:N) {
    
    # Horizontal (x,y)
    if (!is.na(hd[i])) {
      var_xy_i <- var0_xy + var1_xy * hd[i]
      Hmat[i,1] <- var_xy_i
      Hmat[i,2] <- var_xy_i
    } else {
      # No GPS → huge variance
      Hmat[i,1] <- 1e6
      Hmat[i,2] <- 1e6
    }
    
    # Depth
    Hmat[i,3] <- sd_depth^2
  }
  
  Hmat
}






neg_loglikelihood_LinearError2 <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])   # horizontal process noise SD
  sigma2 <- exp(params["sigma2"])   # vertical process noise SD
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  delta <- rep(dt, nrow(data_aug))
  
  # 2-sigma process noise
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  # Build H matrix (linear HDOP)
  var0_xy  <- 0
  var1_xy  <- 0.1
  sd_depth <- 10
  
  Hmat <- build_Hmat_LinearError2(data_aug, var0_xy, var1_xy, sd_depth)
  
  # Initial state
  a <- c(
    ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
    ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
    ifelse(is.na(y[1,3]), 0, y[1,3]), 0
  )
  P <- diag(6) * 1e2
  
  # 2-sigma filter
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

















