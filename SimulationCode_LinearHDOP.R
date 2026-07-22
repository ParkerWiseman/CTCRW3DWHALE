

# Parameter values are not re-obtained.

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