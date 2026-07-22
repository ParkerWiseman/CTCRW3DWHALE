
library(MASS)
library(dplyr)
library(ggplot2)
set.seed(123)

source("matrices.R")

###############################################
# 1. TRUE PARAMETERS
###############################################
beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

s1 <- sigma1_true^2
s2 <- sigma1_true^2
s3 <- sigma2_true^2

###############################################
# 2. TIME GRID
###############################################
N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

###############################################
# 3. LATENT CTCRW SIMULATION
###############################################
X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)



Tmat <- makeT_R(beta1_true, beta2_true, dt)
Qmat <- makeQ_R(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

###############################################
# 4. HDOP SIMULATION
###############################################
hdop <- runif(N, 1, 25)

###############################################
# 5. FILTER + LIKELIHOOD FUNCTIONS
###############################################
makeT_R <- function(b1, b2, delta) {
  e1 <- exp(-b1 * delta)
  e2 <- exp(-b2 * delta)
  T <- matrix(0, 6, 6)
  T[1,1] <- 1; T[2,2] <- e1; T[1,2] <- (1 - e1)/b1
  T[3,3] <- 1; T[4,4] <- e1; T[3,4] <- (1 - e1)/b1
  T[5,5] <- 1; T[6,6] <- e2; T[5,6] <- (1 - e2)/b2
  T
}

makeQ_R <- function(b1, b2, s1, s2, s3, delta) {
  e1  <- exp(-b1 * delta)
  e2  <- exp(-2*b1*delta)
  q11_1 <- (delta - 2*(1-e1)/b1 + (1-e2)/(2*b1)) / b1^2
  q13_1 <- ((1-e1) - (1-e2)/2) / b1^2
  q33_1 <- (1-e2)/(2*b1)
  
  e1v  <- exp(-b2 * delta)
  e2v  <- exp(-2*b2*delta)
  q11_2 <- (delta - 2*(1-e1v)/b2 + (1-e2v)/(2*b2)) / b2^2
  q13_2 <- ((1-e1v) - (1-e2v)/2) / b2^2
  q33_2 <- (1-e2v)/(2*b2)
  
  Q <- matrix(0, 6, 6)
  Q[1,1] <- s1*q11_1; Q[2,2] <- s1*q33_1; Q[1,2] <- Q[2,1] <- s1*q13_1
  Q[3,3] <- s2*q11_1; Q[4,4] <- s2*q33_1; Q[3,4] <- Q[4,3] <- s2*q13_1
  Q[5,5] <- s3*q11_2; Q[6,6] <- s3*q33_2; Q[5,6] <- Q[6,5] <- s3*q13_2
  Q
}

CTCRW_filter <- function(y, Hmat, beta1_vec, beta2_vec,
                         s1, s2, s3, delta, a, P) {
  
  N <- nrow(y)
  
  Z <- matrix(0, 3, 6)
  Z[1,1] <- 1
  Z[2,3] <- 1
  Z[3,5] <- 1
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT_R(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ_R(beta1_vec[i], beta2_vec[i], s1, s2, s3, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    v <- y[i,] - Z %*% aest
    Fmat <- Z %*% Pest %*% t(Z) + diag(Hmat[i,], 3)
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z) + Qmat
  }
  
  ll
}

neg_loglikelihood <- function(params, y, Hmat) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  delta <- rep(dt, nrow(y))
  
  s1 <- sigma1^2
  s2 <- sigma1^2
  s3 <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(y))
  beta2_vec <- rep(beta2, nrow(y))
  
  a <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
  P <- diag(6) * 1e2
  
  ll <- CTCRW_filter(
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
  
  -ll
}

###############################################
# 6. HDOP SENSITIVITY EXPERIMENT
###############################################
hdop_scales <- seq(0, 2, length.out = 20)   # b coefficient from 0 → 2

results <- data.frame(
  b = hdop_scales,
  beta1 = NA,
  beta2 = NA,
  sigma1 = NA,
  sigma2 = NA
)

for (k in seq_along(hdop_scales)) {
  
  b <- hdop_scales[k]
  a <- 5^2
  
  sd_xy <- sqrt(a + b * hdop^2)
  obs_x <- X[,1] + rnorm(N, 0, sd_xy)
  obs_y <- X[,3] + rnorm(N, 0, sd_xy)
  obs_depth <- X[,5] + rnorm(N, 0, 5)
  
  y <- cbind(obs_x, obs_y, obs_depth)
  Hmat <- cbind(a + b * hdop^2, a + b * hdop^2, 5^2)
  
  params_start <- c(
    beta1  = log(1),
    beta2  = log(1),
    sigma1 = log(10),
    sigma2 = log(10)
  )
  
  fit <- optim(
    par      = params_start,
    fn       = neg_loglikelihood,
    y        = y,
    Hmat     = Hmat,
    method   = "L-BFGS-B",
    control  = list(maxit = 500)
  )
  
  p <- exp(fit$par)
  
  results[k,2:5] <- p
}

###############################################
# 7. PLOTS OF PARAMETER DRIFT
###############################################
results_long <- tidyr::pivot_longer(
  results,
  cols = c(beta1, beta2, sigma1, sigma2),
  names_to = "parameter",
  values_to = "estimate"
)

ggplot(results_long, aes(x = b, y = estimate, color = parameter)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal(base_size = 16) +
  labs(
    title = "CTCRW Parameter Drift vs HDOP Measurement Error Scale",
    x = "HDOP Quadratic Coefficient (b)",
    y = "Estimated Parameter Value"
  )