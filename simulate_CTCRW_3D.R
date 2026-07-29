


simulate_CTCRW_3D <- function(
    N,
    dt,
    beta1_true,
    beta2_true,
    sigma1_true,
    sigma2_true,
    x0 = 0,
    y0 = 0,
    depth0 = -50
) {
  # Process noise variances
  s_horiz <- sigma1_true^2
  s_vert  <- sigma2_true^2
  
  # Allocate latent state matrix
  X <- matrix(0, nrow = N, ncol = 6)
  X[1,] <- c(x0, 0, y0, 0, depth0, 0)
  
  # Transition + process noise matrices
  Tmat <- makeT(beta1_true, beta2_true, dt)
  Qmat <- makeQ(beta1_true, beta2_true, s_horiz, s_vert, dt)
  
  # Simulate latent CTCRW
  for (i in 2:N) {
    X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
  }
  
  # Build output data frame
  time_vec <- seq(0, by = dt, length.out = N)
  
  sim_data <- data.frame(
    time  = lubridate::ymd_hms("2020-01-01 00:00:00") + time_vec * 86400,
    x     = X[,1],
    y     = X[,3],
    depth = X[,5]
  )
  
  sim_data
}





