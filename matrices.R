


############################################################
# matrices.R
#
# CTCRW transition and process-noise matrices
#
# State:
#   [x, vx, y, vy, depth, vdepth]
#
# beta1 = horizontal velocity mean-reversion parameter
# beta2 = vertical velocity mean-reversion parameter
#
# sigma/process noise parameters are supplied as variances
# through s_horiz and s_vert.
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
  T[1,2] <- (1 - ebd1) / b1
  
  # Y / horizontal velocity
  T[3,3] <- 1
  T[4,4] <- ebd1
  T[3,4] <- (1 - ebd1) / b1
  
  # Depth / vertical velocity
  T[5,5] <- 1
  T[6,6] <- ebd2
  T[5,6] <- (1 - ebd2) / b2
  
  T
}

############################################################
# PROCESS NOISE MATRIX
############################################################

makeQ <- function(
    b1,
    b2,
    s_horiz,
    s_vert,
    delta
) {
  
  ##########################################################
  # Horizontal component
  ##########################################################
  
  ebd1  <- exp(-b1 * delta)
  e2bd1 <- exp(-2 * b1 * delta)
  
  q11_1 <-
    (
      delta -
        2 * (1 - ebd1) / b1 +
        (1 - e2bd1) / (2 * b1)
    ) / b1^2
  
  q13_1 <-
    (
      (1 - ebd1) -
        (1 - e2bd1) / 2
    ) / b1^2
  
  q33_1 <-
    (1 - e2bd1) / (2 * b1)
  
  ##########################################################
  # Vertical component
  ##########################################################
  
  ebd2  <- exp(-b2 * delta)
  e2bd2 <- exp(-2 * b2 * delta)
  
  q11_2 <-
    (
      delta -
        2 * (1 - ebd2) / b2 +
        (1 - e2bd2) / (2 * b2)
    ) / b2^2
  
  q13_2 <-
    (
      (1 - ebd2) -
        (1 - e2bd2) / 2
    ) / b2^2
  
  q33_2 <-
    (1 - e2bd2) / (2 * b2)
  
  ##########################################################
  # Assemble 6 x 6 Q matrix
  ##########################################################
  
  Q <- matrix(0, 6, 6)
  
  # X
  Q[1,1] <- s_horiz * q11_1
  Q[2,2] <- s_horiz * q33_1
  Q[1,2] <- s_horiz * q13_1
  Q[2,1] <- s_horiz * q13_1
  
  # Y
  Q[3,3] <- s_horiz * q11_1
  Q[4,4] <- s_horiz * q33_1
  Q[3,4] <- s_horiz * q13_1
  Q[4,3] <- s_horiz * q13_1
  
  # Depth
  Q[5,5] <- s_vert * q11_2
  Q[6,6] <- s_vert * q33_2
  Q[5,6] <- s_vert * q13_2
  Q[6,5] <- s_vert * q13_2
  
  Q
}

############################################################
# OPTIONAL ALIASES
#
# These allow old code using makeT_R / makeQ_R to continue
# working.
############################################################

makeT_R <- function(b1, b2, delta) {
  
  makeT(
    b1 = b1,
    b2 = b2,
    delta = delta
  )
}

makeQ_R <- function(b1, b2, s1, s2, s3, delta) {
  
  Q <- matrix(0, 6, 6)
  
  ##########################################################
  # Horizontal
  ##########################################################
  
  ebd1  <- exp(-b1 * delta)
  e2bd1 <- exp(-2 * b1 * delta)
  
  q11_1 <-
    (
      delta -
        2 * (1 - ebd1) / b1 +
        (1 - e2bd1) / (2 * b1)
    ) / b1^2
  
  q13_1 <-
    (
      (1 - ebd1) -
        (1 - e2bd1) / 2
    ) / b1^2
  
  q33_1 <-
    (1 - e2bd1) / (2 * b1)
  
  Q[1,1] <- s1 * q11_1
  Q[2,2] <- s1 * q33_1
  Q[1,2] <- s1 * q13_1
  Q[2,1] <- s1 * q13_1
  
  Q[3,3] <- s2 * q11_1
  Q[4,4] <- s2 * q33_1
  Q[3,4] <- s2 * q13_1
  Q[4,3] <- s2 * q13_1
  
  ##########################################################
  # Vertical
  ##########################################################
  
  ebd2  <- exp(-b2 * delta)
  e2bd2 <- exp(-2 * b2 * delta)
  
  q11_2 <-
    (
      delta -
        2 * (1 - ebd2) / b2 +
        (1 - e2bd2) / (2 * b2)
    ) / b2^2
  
  q13_2 <-
    (
      (1 - ebd2) -
        (1 - e2bd2) / 2
    ) / b2^2
  
  q33_2 <-
    (1 - e2bd2) / (2 * b2)
  
  Q[5,5] <- s3 * q11_2
  Q[6,6] <- s3 * q33_2
  Q[5,6] <- s3 * q13_2
  Q[6,5] <- s3 * q13_2
  
  Q
}














