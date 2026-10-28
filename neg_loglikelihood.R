


############################################################
# NEGATIVE LOG-LIKELIHOOD (small measurement error)
############################################################


source("utility.R")


neg_loglikelihood_noerror <- function(
    params,
    data_aug
) {
  
  beta_horiz  <- exp(params["beta_horiz"])
  beta_vert   <- exp(params["beta_vert"])
  sigma_horiz <- exp(params["sigma_horiz"])
  sigma_vert  <- exp(params["sigma_vert"])
  
  y <- as.matrix(data_aug[, c("x", "y", "depth")])
  N <- nrow(data_aug)
  
  # Using the numeric time column Time_sec directly (and removing two unnecessary lines)
  delta <- c(diff(data_aug$Time_sec),0)
  
  var_horiz_vec <- rep(sigma_horiz^2, N)
  var_vert_vec  <- rep(sigma_vert^2,  N)
  
  beta_horiz_vec <- rep(beta_horiz, N)
  beta_vert_vec  <- rep(beta_vert,  N)
  
  Hmat <- matrix(0, N, 3)
  Hmat[,1] <- 5^2   # x
  Hmat[,2] <- 5^2   # y
  Hmat[,3] <- 10^2  # depth
  
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  #P <- diag(6) * 1e6
  P <- diag(c(100, 10, 100, 10, 100, 10))
  
  filt <- CTCRW_filter(
    y              = y,
    Hmat           = Hmat,
    beta_horiz_vec = beta_horiz_vec,
    beta_vert_vec  = beta_vert_vec,
    var_horiz_vec  = var_horiz_vec,
    var_vert_vec   = var_vert_vec,
    delta          = delta,
    a              = a,
    P              = P
  )
  
  -filt$ll
}









