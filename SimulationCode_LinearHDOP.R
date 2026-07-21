

# Parameter values are not re-obtained.

library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

# TRUE PARAMETERS
beta1_true  <- 2.8     # horizontal OU damping
beta2_true  <- 0.8     # vertical OU damping
sigma1_true <- 30      # horizontal process noise SD
sigma2_true <- 10      # vertical process noise SD

s1 <- sigma1_true^2
s2 <- sigma1_true^2
s3 <- sigma2_true^2

# TIME GRID
N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

# LATENT STATE: (x, vx, y, vy, depth, vdepth)
X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)

makeT <- function(b1, b2, dt) {
  T <- matrix(0, 6, 6)
  e1 <- exp(-b1 * dt)
  e2 <- exp(-b2 * dt)
  T[1,1] <- 1; T[2,2] <- e1; T[1,2] <- (1 - e1)/b1
  T[3,3] <- 1; T[4,4] <- e1; T[3,4] <- (1 - e1)/b1
  T[5,5] <- 1; T[6,6] <- e2; T[5,6] <- (1 - e2)/b2
  T
}

makeQ <- function(b1, b2, s1, s2, s3, dt) {
  e1  <- exp(-b1 * dt)
  e2  <- exp(-2*b1*dt)
  q11_1 <- (dt - 2*(1-e1)/b1 + (1-e2)/(2*b1)) / b1^2
  q13_1 <- ((1-e1) - (1-e2)/2) / b1^2
  q33_1 <- (1-e2)/(2*b1)
  
  e1v  <- exp(-b2 * dt)
  e2v  <- exp(-2*b2*dt)
  q11_2 <- (dt - 2*(1-e1v)/b2 + (1-e2v)/(2*b2)) / b2^2
  q13_2 <- ((1-e1v) - (1-e2v)/2) / b2^2
  q33_2 <- (1-e2v)/(2*b2)
  
  Q <- matrix(0, 6, 6)
  Q[1,1] <- s1*q11_1; Q[2,2] <- s1*q33_1; Q[1,2] <- Q[2,1] <- s1*q13_1
  Q[3,3] <- s2*q11_1; Q[4,4] <- s2*q33_1; Q[3,4] <- Q[4,3] <- s2*q13_1
  Q[5,5] <- s3*q11_2; Q[6,6] <- s3*q33_2; Q[5,6] <- Q[6,5] <- s3*q13_2
  Q
}

# SIMULATE LATENT PROCESS
Tmat <- makeT(beta1_true, beta2_true, dt)
Qmat <- makeQ(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

# OBSERVATIONS WITH LINEAR HDOP MEASUREMENT ERROR
obs_type <- rbinom(N, 1, 0.7)

obs_x     <- rep(NA, N)
obs_y     <- rep(NA, N)
obs_depth <- rep(NA, N)
hdop      <- rep(NA, N)

# LINEAR HDOP MODEL:
# Var(x) = var0_xy + var1_xy * HDOP
var0_xy <- 0
var1_xy <- 0.1
# Changing these two values changes how well the re-obtained values match the chosen values.
sd_depth  <- 10

sd_xy_fun <- function(h) {
  sqrt(var0_xy + var1_xy * h)   # LINEAR HDOP
}

for (i in 1:N) {
  if (obs_type[i] == 1) {
    obs_depth[i] <- X[i,5] + rnorm(1, 0, sd_depth)
    hdop[i]      <- NA
  } else {
    hdop[i]  <- runif(1, 1, 25)
    sd_xy_i  <- sd_xy_fun(hdop[i])
    obs_x[i] <- X[i,1] + rnorm(1, 0, sd_xy_i)
    obs_y[i] <- X[i,3] + rnorm(1, 0, sd_xy_i)
  }
}

sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth,
  hdop  = hdop
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    hasObs = !is.na(x) | !is.na(y) | !is.na(depth),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

# FILTER MATRICES
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

# BUILD Hmat USING LINEAR HDOP
build_Hmat <- function(data_aug, var0_xy, var1_xy, sd_depth) {
  y  <- as.matrix(data_aug[, c("x","y","depth")])
  hd <- data_aug$hdop
  
  Hmat <- matrix(NA_real_, nrow(data_aug), 3)
  
  for (i in 1:nrow(data_aug)) {
    if (!is.na(hd[i])) {
      var_xy_i <- var0_xy + var1_xy * hd[i]   # LINEAR HDOP
      if (!is.na(y[i,1])) Hmat[i,1] <- var_xy_i
      if (!is.na(y[i,2])) Hmat[i,2] <- var_xy_i
    }
    if (!is.na(y[i,3])) {
      Hmat[i,3] <- sd_depth^2
    }
  }
  
  Hmat
}

# KALMAN FILTER
CTCRW_filter <- function(y, Hmat, beta1_vec, beta2_vec,
                         s1, s2, s3, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT_R(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ_R(beta1_vec[i], beta2_vec[i], s1, s2, s3, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    obs_idx <- which(!is.na(y[i, ]))
    
    if (length(obs_idx) == 0) {
      aest <- a_pred
      Pest <- P_pred
      a_f[i,]  <- as.numeric(aest)
      P_f[[i]] <- Pest
      next
    }
    
    Z_i <- Z_full[obs_idx, , drop = FALSE]
    H_i <- diag(Hmat[i, obs_idx], length(obs_idx))
    
    v <- y[i, obs_idx] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    detF <- det(Fmat)
    if (!is.finite(detF) || detF <= 0) {
      aest <- a_pred
      Pest <- P_pred
      a_f[i,]  <- as.numeric(aest)
      P_f[[i]] <- Pest
      next
    }
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(detF) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,]  <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}

# SMOOTHER
CTCRW_smoother <- function(filter_out, beta1_vec, beta2_vec,
                           s1, s2, s3, delta) {
  
  N <- nrow(filter_out$a_f)
  a_s <- filter_out$a_f
  P_s <- filter_out$P_f
  
  for (i in (N-1):1) {
    
    Tmat <- makeT_R(beta1_vec[i+1], beta2_vec[i+1], delta[i+1])
    
    a_f_i   <- filter_out$a_f[i,]
    a_p_ip1 <- filter_out$a_p[i+1,]
    P_f_i   <- filter_out$P_f[[i]]
    P_p_ip1 <- filter_out$P_p[[i+1]]
    
    J <- tryCatch(P_f_i %*% t(Tmat) %*% solve(P_p_ip1), error=function(e) NULL)
    if (is.null(J)) next
    
    a_s[i,]  <- a_f_i + J %*% (a_s[i+1,] - a_p_ip1)
    P_s[[i]] <- P_f_i + J %*% (P_s[[i+1]] - P_p_ip1) %*% t(J)
  }
  
  list(a_s = a_s, P_s = P_s)
}

# NEGATIVE LOG-LIKELIHOOD
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

params_start <- c(
  beta1  = log(0.5),
  beta2  = log(0.5),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
p_hat