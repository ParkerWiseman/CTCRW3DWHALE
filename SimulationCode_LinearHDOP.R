

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









params_start <- c(
  beta1  = log(0.5),
  beta2  = log(0.5),
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




####################################
####################################





# (using s_horiz and s_vert)






library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

source("matrices.R")          # must contain makeT() and makeQ()
source("CTCRW_filter.R")      # must contain CTCRW_filter1()
source("CTCRW_smoother.R")    # must contain CTCRW_smoother1()
source("neg_loglikelihood.R") # we will define new functions below

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30   # horizontal
sigma2_true <- 10   # vertical

s_horiz_true <- sigma1_true^2
s_vert_true  <- sigma2_true^2

###############################################
# SIMULATION (2-sigma CTCRW)
###############################################

N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)

Tmat <- makeT(beta1_true, beta2_true, dt)
Qmat <- makeQ(beta1_true, beta2_true, s_horiz_true, s_vert_true, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

###############################################
# LINEAR HDOP MEASUREMENT ERROR MODEL
###############################################

obs_type <- rbinom(N, 1, 0.7)

obs_x     <- rep(NA, N)
obs_y     <- rep(NA, N)
obs_depth <- rep(NA, N)
hdop      <- rep(NA, N)

var0_xy <- 0
var1_xy <- 0.1
sd_depth <- 10

sd_xy_fun <- function(h) sqrt(var0_xy + var1_xy * h)

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
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])





###############################################
# OPTIMIZATION
###############################################

params_start <- c(
  beta1  = log(0.5),
  beta2  = log(0.5),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_LinearError2,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

exp(fit$par)








#########################################




# IGNORE ALL ABOVE CODE
# GOOD CODE
# uses s_horiz and s_vert, and simulation function





library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

source("matrices.R")          # makeT(), makeQ()
source("CTCRW_filter.R")      # CTCRW_filter1()
source("CTCRW_smoother.R")    # CTCRW_smoother1()
source("neg_loglikelihood.R") # neg_loglikelihood_LinearError2()
source("simulate_CTCRW_3D.R") # NEW simulation function

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

###############################################
# SIMULATION USING THE NEW FUNCTION
###############################################

N  <- 5000
dt <- 15 / (24*60)

sim_latent <- simulate_CTCRW_3D(
  N           = N,
  dt          = dt,
  beta1_true  = beta1_true,
  beta2_true  = beta2_true,
  sigma1_true = sigma1_true,
  sigma2_true = sigma2_true
)

###############################################
# LINEAR HDOP MEASUREMENT ERROR MODEL
###############################################

obs_type <- rbinom(N, 1, 0.7)

obs_x     <- rep(NA, N)
obs_y     <- rep(NA, N)
obs_depth <- rep(NA, N)
hdop      <- rep(NA, N)

var0_xy <- 0
var1_xy <- 0.1
sd_depth <- 10

sd_xy_fun <- function(h) sqrt(var0_xy + var1_xy * h)

for (i in 1:N) {
  if (obs_type[i] == 1) {
    obs_depth[i] <- sim_latent$depth[i] + rnorm(1, 0, sd_depth)
    hdop[i]      <- NA
  } else {
    hdop[i]  <- runif(1, 1, 25)
    sd_xy_i  <- sd_xy_fun(hdop[i])
    obs_x[i] <- sim_latent$x[i] + rnorm(1, 0, sd_xy_i)
    obs_y[i] <- sim_latent$y[i] + rnorm(1, 0, sd_xy_i)
  }
}

sim_data <- data.frame(
  time  = sim_latent$time,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth,
  hdop  = hdop
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

###############################################
# OPTIMIZATION
###############################################

params_start <- c(
  beta1  = log(0.5),
  beta2  = log(0.5),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_LinearError2,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

exp(fit$par)






######################################

# August 12, 2026




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
source("simulate_CTCRW_3D.R") 

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

###############################################
# SIMULATION USING THE NEW FUNCTION
###############################################

N  <- 5000
dt <- 15 / (24*60)

sim_latent <- simulate_CTCRW_3D(
  N           = N,
  dt          = dt,
  beta1_true  = beta1_true,
  beta2_true  = beta2_true,
  sigma1_true = sigma1_true,
  sigma2_true = sigma2_true
)

###############################################
# LINEAR HDOP MEASUREMENT ERROR MODEL
###############################################

obs_type <- rbinom(N, 1, 0.7)

obs_x     <- rep(NA, N)
obs_y     <- rep(NA, N)
obs_depth <- rep(NA, N)
hdop      <- rep(NA, N)

var0_xy  <- 0
var1_xy  <- 0.1
sd_depth <- 10

sd_xy_fun <- function(h) sqrt(var0_xy + var1_xy * h)

for (i in 1:N) {
  if (obs_type[i] == 1) {
    obs_depth[i] <- sim_latent$depth[i] + rnorm(1, 0, sd_depth)
    hdop[i]      <- NA
  } else {
    hdop[i]  <- runif(1, 1, 25)
    sd_xy_i  <- sd_xy_fun(hdop[i])
    obs_x[i] <- sim_latent$x[i] + rnorm(1, 0, sd_xy_i)
    obs_y[i] <- sim_latent$y[i] + rnorm(1, 0, sd_xy_i)
  }
}

sim_data <- data.frame(
  time  = sim_latent$time,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth,
  hdop  = hdop
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

###############################################
# OPTIMIZATION USING UNIFIED LIKELIHOOD
###############################################

params_start <- c(
  beta1  = log(0.5),
  beta2  = log(0.5),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood(p, aug, error_model = "linearerror"),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

exp(fit$par)











###################################

# histograms from multiple simulations








library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")
source("simulate_CTCRW_3D.R")

set.seed(123)

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

###############################################
# MONTE CARLO SETTINGS
###############################################

n_sims <- 50
N      <- 5000
dt     <- 15 / (24*60)

# storage for parameter estimates
estimates <- data.frame(
  beta1  = numeric(n_sims),
  beta2  = numeric(n_sims),
  sigma1 = numeric(n_sims),
  sigma2 = numeric(n_sims)
)

###############################################
# MONTE CARLO LOOP
###############################################

for (s in 1:n_sims) {
  
  # --- Simulate latent CTCRW ---
  sim_latent <- simulate_CTCRW_3D(
    N           = N,
    dt          = dt,
    beta1_true  = beta1_true,
    beta2_true  = beta2_true,
    sigma1_true = sigma1_true,
    sigma2_true = sigma2_true
  )
  
  # --- Linear HDOP measurement error model ---
  obs_type <- rbinom(N, 1, 0.7)
  
  obs_x     <- rep(NA, N)
  obs_y     <- rep(NA, N)
  obs_depth <- rep(NA, N)
  hdop      <- rep(NA, N)
  
  var0_xy  <- 0
  var1_xy  <- 0.1
  sd_depth <- 10
  
  sd_xy_fun <- function(h) sqrt(var0_xy + var1_xy * h)
  
  for (i in 1:N) {
    if (obs_type[i] == 1) {
      obs_depth[i] <- sim_latent$depth[i] + rnorm(1, 0, sd_depth)
      hdop[i]      <- NA
    } else {
      hdop[i]  <- runif(1, 1, 25)
      sd_xy_i  <- sd_xy_fun(hdop[i])
      obs_x[i] <- sim_latent$x[i] + rnorm(1, 0, sd_xy_i)
      obs_y[i] <- sim_latent$y[i] + rnorm(1, 0, sd_xy_i)
    }
  }
  
  sim_data <- data.frame(
    time  = sim_latent$time,
    x     = obs_x,
    y     = obs_y,
    depth = obs_depth,
    hdop  = hdop
  )
  
  aug <- sim_data %>%
    mutate(
      Time = as.numeric(difftime(time, min(time), units = "days")),
      orig_index = seq_len(n())
    )
  
  # --- Optimization ---
  params_start <- c(
    beta1  = log(0.5),
    beta2  = log(0.5),
    sigma1 = log(10),
    sigma2 = log(10)
  )
  
  fit <- optim(
    par      = params_start,
    fn       = function(p) neg_loglikelihood(p, aug, error_model = "linearerror"),
    method   = "L-BFGS-B",
    control  = list(trace = 0, maxit = 1000)
  )
  
  p_hat <- exp(fit$par)
  
  # store estimates
  estimates[s, ] <- p_hat
}

###############################################
# HISTOGRAMS OF PARAMETER ESTIMATES
###############################################

par(mfrow=c(2,2))

hist(estimates$beta1, main="β1 estimates", xlab="beta1", col="skyblue")
abline(v=beta1_true, col="red", lwd=2)

hist(estimates$beta2, main="β2 estimates", xlab="beta2", col="skyblue")
abline(v=beta2_true, col="red", lwd=2)

hist(estimates$sigma1, main="σ1 estimates", xlab="sigma1", col="skyblue")
abline(v=sigma1_true, col="red", lwd=2)

hist(estimates$sigma2, main="σ2 estimates", xlab="sigma2", col="skyblue")
abline(v=sigma2_true, col="red", lwd=2)

par(mfrow=c(1,1))

###############################################
# SUMMARY STATISTICS
###############################################

summary(estimates)




















