


############################################################
# simulate_CTCRW_3D.R
#
# Simulate a 3D CTCRW trajectory
#
# IMPORTANT:
#
# beta1 and beta2 are in 1/second.
# sigma1 and sigma2 are in movement units / sqrt(second)
#
# Therefore the CTCRW transition matrices must receive
# dt in SECONDS.
#
# The output time column is converted to POSIXct separately.
############################################################

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
  
  
  ##########################################################
  # IMPORTANT:
  #
  # dt is expected to be in SECONDS.
  #
  # Example:
  #
  # dt = 60
  #
  # means one state every 60 seconds.
  ##########################################################
  
  
  ##########################################################
  # PROCESS NOISE VARIANCES
  ##########################################################
  
  s_horiz <- sigma1_true^2
  
  s_vert <- sigma2_true^2
  
  
  ##########################################################
  # ALLOCATE LATENT STATE MATRIX
  #
  # State:
  #
  # 1 = x
  # 2 = vx
  # 3 = y
  # 4 = vy
  # 5 = depth
  # 6 = vertical velocity
  ##########################################################
  
  X <- matrix(
    0,
    nrow = N,
    ncol = 6
  )
  
  
  ##########################################################
  # INITIAL STATE
  ##########################################################
  
  X[1,] <- c(
    x0,
    0,
    y0,
    0,
    depth0,
    0
  )
  
  
  ##########################################################
  # TRANSITION MATRIX
  #
  # IMPORTANT:
  #
  # dt is now in SECONDS.
  ##########################################################
  
  Tmat <- makeT(
    beta1_true,
    beta2_true,
    dt
  )
  
  
  ##########################################################
  # PROCESS NOISE MATRIX
  #
  # IMPORTANT:
  #
  # dt is now in SECONDS.
  ##########################################################
  
  Qmat <- makeQ(
    beta1_true,
    beta2_true,
    s_horiz,
    s_vert,
    dt
  )
  
  
  ##########################################################
  # SIMULATE CTCRW
  ##########################################################
  
  for (i in 2:N) {
    
    X[i,] <-
      Tmat %*% X[i-1,] +
      MASS::mvrnorm(
        1,
        rep(0, 6),
        Qmat
      )
  }
  
  
  ##########################################################
  # CREATE TIME VECTOR
  #
  # POSIXct requires seconds.
  #
  # We therefore use:
  #
  # (0, dt, 2*dt, ...) seconds
  #
  # directly.
  ##########################################################
  
  time_seconds <-
    seq(
      0,
      by = dt,
      length.out = N
    )
  
  
  ##########################################################
  # CREATE OUTPUT DATA FRAME
  ##########################################################
  
  sim_data <- data.frame(
    
    time =
      as.POSIXct(
        "2020-01-01 00:00:00",
        tz = "UTC"
      ) +
      time_seconds,
    
    x = X[,1],
    
    y = X[,3],
    
    depth = X[,5]
  )
  
  
  ##########################################################
  # RETURN
  ##########################################################
  
  sim_data
}











