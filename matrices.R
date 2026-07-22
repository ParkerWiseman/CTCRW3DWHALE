

makeT_R <- function(b1, b2, delta) {
  ebd1 <- exp(-b1 * delta)
  ebd2 <- exp(-b2 * delta)
  T <- matrix(0, 6, 6)
  T[1,1] <- 1; T[2,2] <- ebd1; T[1,2] <- (1 - ebd1)/b1
  T[3,3] <- 1; T[4,4] <- ebd1; T[3,4] <- (1 - ebd1)/b1
  T[5,5] <- 1; T[6,6] <- ebd2; T[5,6] <- (1 - ebd2)/b2
  T
}





makeQ_R <- function(b1, b2, s1, s2, s3, delta) {
  ebd1  <- exp(-b1 * delta)
  e2bd1 <- exp(-2*b1*delta)
  q11_1 <- (delta - 2*(1-ebd1)/b1 + (1-e2bd1)/(2*b1)) / b1^2
  q13_1 <- ((1-ebd1) - (1-e2bd1)/2) / b1^2
  q33_1 <- (1-e2bd1)/(2*b1)
  
  ebd2  <- exp(-b2 * delta)
  e2bd2 <- exp(-2*b2*delta)
  q11_2 <- (delta - 2*(1-ebd2)/b2 + (1-e2bd2)/(2*b2)) / b2^2
  q13_2 <- ((1-ebd2) - (1-e2bd2)/2) / b2^2
  q33_2 <- (1-e2bd2)/(2*b2)
  
  Q <- matrix(0, 6, 6)
  Q[1,1] <- s1*q11_1; Q[2,2] <- s1*q33_1; Q[1,2] <- Q[2,1] <- s1*q13_1
  Q[3,3] <- s2*q11_1; Q[4,4] <- s2*q33_1; Q[3,4] <- Q[4,3] <- s2*q13_1
  Q[5,5] <- s3*q11_2; Q[6,6] <- s3*q33_2; Q[5,6] <- Q[6,5] <- s3*q13_2
  Q
}















