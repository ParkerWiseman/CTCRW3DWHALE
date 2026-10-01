


############################################################
# matrices.R
#
# CTCRW transition and process-noise matrices
#
# State:
#   [x, vx, y, vy, depth, vdepth]
#
# beta_horiz = horizontal velocity mean-reversion parameter
# beta_vert = vertical velocity mean-reversion parameter
#
# sigma/process noise parameters are supplied as variances
# through var_horiz_vec and var_vert_vec.
############################################################


############################################################
# TRANSITION MATRIX
############################################################

makeT <- function(b1, b2, delta) {
  
  ebd1 <- exp(-b1 * delta)
  ebd2 <- exp(-b2 * delta)
  
  T <- matrix(0, 6, 6)
  
  # X / horizontal velocity
  T[1,1] <- 1
  T[2,2] <- ebd1
  T[1,2] <- (1-ebd1)/b1
  
  # Y / horizontal velocity
  T[3,3] <- 1
  T[4,4] <- ebd1
  T[3,4] <- (1-ebd1)/b1
  
  # Depth / vertical velocity
  T[5,5] <- 1
  T[6,6] <- ebd2
  T[5,6] <- (1-ebd2)/b2
  
  T
}

############################################################
# PROCESS NOISE MATRIX
############################################################

makeQ <- function(
    b1,
    b2,
    var_horiz,
    var_vert,
    delta
) {
  
  ##########################################################
  # Horizontal component
  ##########################################################
  
  ebd1  <- exp(-b1*delta)
  e2bd1 <- exp(-2*b1*delta)
  q11_1 <- (delta - 2*(1-ebd1)/b1 + (1-e2bd1)/(2*b1)) / b1^2
  q13_1 <- ((1-ebd1) - (1-e2bd1)/2) / b1^2
  q33_1 <- (1-e2bd1)/(2*b1)
  
  ##########################################################
  # Vertical component
  ##########################################################
  
  ebd2  <- exp(-b2*delta)
  e2bd2 <- exp(-2*b2*delta)
  q11_2 <- (delta - 2*(1-ebd2)/b2 + (1-e2bd2)/(2*b2)) / b2^2
  q13_2 <- ((1-ebd2) - (1-e2bd2)/2) / b2^2
  q33_2 <- (1-e2bd2)/(2*b2)
  
  ##########################################################
  # Assemble 6 x 6 Q matrix
  ##########################################################
  
  Q <- matrix(0, 6, 6)
  
  # X
  Q[1,1] <- var_horiz * q11_1
  Q[2,2] <- var_horiz * q33_1
  Q[1,2] <- var_horiz * q13_1
  Q[2,1] <- var_horiz * q13_1
  
  # Y
  Q[3,3] <- var_horiz * q11_1
  Q[4,4] <- var_horiz * q33_1
  Q[3,4] <- var_horiz * q13_1
  Q[4,3] <- var_horiz * q13_1
  
  # Depth
  Q[5,5] <- var_vert * q11_2
  Q[6,6] <- var_vert * q33_2
  Q[5,6] <- var_vert * q13_2
  Q[6,5] <- var_vert * q13_2
  
  Q
}






