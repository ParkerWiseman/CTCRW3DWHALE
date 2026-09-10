


##########################

# September 6, 2026

# using Michelot & Blackwell formulas AND USING HDOP=0 FOR ALL OBSERVATIONS

############################################################
# neg_loglikelihood.R — CTCRW likelihood with delta in seconds
############################################################

build_Hmat_LinearError <- function(data_aug, var0_xy, var1_xy, sd_depth) {
  
  hd <- data_aug$hdop
  N  <- nrow(data_aug)
  
  Hmat <- matrix(NA_real_, N, 3)
  
  for (i in 1:N) {
    
    h_i <- hd[i]
    if (is.na(h_i)) h_i <- 10
    
    var_xy_i <- var0_xy + var1_xy * h_i
    var_xy_i <- max(var_xy_i, 1e-3)
    
    Hmat[i,1] <- var_xy_i
    Hmat[i,2] <- var_xy_i
    Hmat[i,3] <- sd_depth^2
  }
  
  Hmat
}

neg_loglikelihood <- function(params, data_aug,
                              error_model = c("linearerror")) {
  
  error_model <- match.arg(error_model)
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  #### ----------------------------------------------------
  #### DELTA IN SECONDS (minimum = 1 second)
  #### ----------------------------------------------------
  
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  
  delta_raw <- c(NA, diff(time_sec))
  delta_raw[1] <- delta_raw[2]
  
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  #### ----------------------------------------------------
  #### PROCESS NOISE
  #### ----------------------------------------------------
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### ----------------------------------------------------
  #### H MATRIX (your chosen values)
  #### ----------------------------------------------------
  
  var0_xy  <- 100
  var1_xy  <- 2
  sd_depth <- 25
  
  Hmat <- build_Hmat_LinearError(data_aug, var0_xy, var1_xy, sd_depth)
  
  #### ----------------------------------------------------
  #### INITIAL STATE (first non-missing)
  #### ----------------------------------------------------
  
  get_first_non_missing <- function(col) {
    idx <- which(!is.na(col))[1]
    if (is.na(idx)) 0 else col[idx]
  }
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  #### ----------------------------------------------------
  #### INITIAL COVARIANCE
  #### ----------------------------------------------------
  
  P <- diag(6) * 1e6
  
  #### ----------------------------------------------------
  #### FILTER
  #### ----------------------------------------------------
  
  filt <- CTCRW_filter1(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz   = s_horiz,
    s_vert    = s_vert,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  -filt$ll
}








############################################################
# CTCRW 3D SPERM WHALE MODEL — delta in seconds + physics-based initials
############################################################

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

############################################################
# LOAD DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")
whale$hdop <- 0
whale$time <- ymd_hms(whale$time)

aug <- whale %>%
  mutate(
    Time_sec = as.numeric(difftime(time, min(time), units = "secs")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

############################################################
# INITIAL VALUES FROM MICHELOT & BLACKWELL (2020)
############################################################

## Autocorrelation times (hours)
tau_horiz_hours <- 8
tau_vert_hours  <- 3

## Convert to seconds
tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours  * 3600

## Beta
beta1_start <- 1 / tau_horiz_sec
beta2_start <- 1 / tau_vert_sec

## RMS speeds (m/s)
v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

## Sigma
sigma1_start <- v_rms_horiz * sqrt(2 * beta1_start)
sigma2_start <- v_rms_vert  * sqrt(2 * beta2_start)

## Log-scale parameters for optim
params_start <- c(
  beta1  = log(beta1_start),
  beta2  = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

print(params_start)

############################################################
# FIT MODEL
############################################################

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood(p, aug, error_model = "linearerror"),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

s_horiz_hat <- sigma1_hat^2
s_vert_hat  <- sigma2_hat^2

############################################################
# BUILD H MATRIX
############################################################

var0_xy  <- 100
var1_xy  <- 2
sd_depth <- 25

Hmat <- build_Hmat_LinearError(
  aug,
  var0_xy = var0_xy,
  var1_xy = var1_xy,
  sd_depth = sd_depth
)

############################################################
# DELTA (seconds)
############################################################

delta_raw <- c(NA, diff(aug$Time_sec))
delta_raw[1] <- delta_raw[2]

delta <- pmax(delta_raw, 1)
delta[!is.finite(delta)] <- 1

############################################################
# INITIAL STATE
############################################################

get_first_non_missing <- function(col) {
  idx <- which(!is.na(col))[1]
  if (is.na(idx)) 0 else col[idx]
}

a0 <- c(
  get_first_non_missing(y[,1]), 0,
  get_first_non_missing(y[,2]), 0,
  get_first_non_missing(y[,3]), 0
)

P0 <- diag(6) * 1e6

############################################################
# FILTER
############################################################

filt <- CTCRW_filter1(
  y         = y,
  Hmat      = Hmat,
  beta1_vec = rep(beta1_hat, N),
  beta2_vec = rep(beta2_hat, N),
  s_horiz   = s_horiz_hat,
  s_vert    = s_vert_hat,
  delta     = delta,
  a         = a0,
  P         = P0
)

############################################################
# SMOOTHER
############################################################

smooth <- CTCRW_smoother1(
  filter_out = filt,
  beta1_vec  = rep(beta1_hat, N),
  beta2_vec  = rep(beta2_hat, N),
  s_horiz    = s_horiz_hat,
  s_vert     = s_vert_hat,
  delta      = delta
)

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time

############################################################
# PLOTS
############################################################

## 1. XY plot
ggplot() +
  geom_point(data = aug, aes(x = x, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = y), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (XY)",
       x = "X", y = "Y") +
  theme_minimal()

## 2. Depth vs Time
ggplot() +
  geom_point(data = aug, aes(x = time, y = depth), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = time, y = depth), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Depth vs Time)",
       x = "Time", y = "Depth") +
  theme_minimal()

## 3. X vs Time
ggplot() +
  geom_point(data = aug, aes(x = time, y = x), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = time, y = x), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (X vs Time)",
       x = "Time", y = "X") +
  theme_minimal()

## 4. Y vs Time
ggplot() +
  geom_point(data = aug, aes(x = time, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = time, y = y), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Y vs Time)",
       x = "Time", y = "Y") +
  theme_minimal()

## 5. 3D plot
plot_ly() %>%
  add_markers(
    data = aug,
    x = ~x, y = ~y, z = ~depth,
    marker = list(color = 'blue', size = 2),
    name = "Observed"
  ) %>%
  add_lines(
    data = smooth_track,
    x = ~x, y = ~y, z = ~depth,
    line = list(color = 'red', width = 4),
    name = "Smoothed"
  ) %>%
  layout(
    title = "3D Whale Track: Observed vs Smoothed",
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Depth")
    )
  )












