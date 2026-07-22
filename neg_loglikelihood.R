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