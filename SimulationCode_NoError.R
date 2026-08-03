# good code
library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")



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



# SIMULATE LATENT PROCESS
Tmat <- makeT_R(beta1_true, beta2_true, dt)
Qmat <- makeQ_R(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

# OBSERVATIONS = TRUE LATENT STATES (NO MEASUREMENT ERROR)
sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = X[,1],
  y     = X[,3],
  depth = X[,5]
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])




params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_NoError,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
p_hat






################
################










# (using just s_horiz and s_vert)


library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")



# TRUE PARAMETERS
beta1_true  <- 2.8     # horizontal OU damping
beta2_true  <- 0.8     # vertical OU damping
sigma1_true <- 30      # horizontal process noise SD
sigma2_true <- 10      # vertical process noise SD


s_horiz <- sigma1_true^2
s_vert  <- sigma2_true^2

# TIME GRID
N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

# LATENT STATE: (x, vx, y, vy, depth, vdepth)
X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)



# SIMULATE LATENT PROCESS
Tmat <- makeT(beta1_true, beta2_true, dt)
Qmat <- makeQ(beta1_true, beta2_true, s_horiz, s_vert, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

# OBSERVATIONS = TRUE LATENT STATES (NO MEASUREMENT ERROR)
sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = X[,1],
  y     = X[,3],
  depth = X[,5]
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])




params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_NoError1,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
p_hat




#########################################

















# GOOD CODE
# uses s_horiz and s_vert, and simulation function



library(MASS)
library(dplyr)
library(lubridate)

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")
source("simulate_CTCRW_3D.R")

# --- Simulation ---
N  <- 5000
dt <- 15 / (24*60)

sim_data <- simulate_CTCRW_3D(
  N          = N,
  dt         = dt,
  beta1_true = 2.8,
  beta2_true = 0.8,
  sigma1_true = 30,
  sigma2_true = 10
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

# --- Optimization ---
params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_NoError1,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

exp(fit$par)












