


############################################################
# SpermWhale_NoMeasurementError_Good.R
#
# 3D CTCRW sperm whale model
#
# Estimates:
#   beta1
#   beta2
#   sigma1
#   sigma2
#
# Then performs:
#   1. Kalman filtering
#   2. Kalman smoothing
#   3. Plotting of observations + smoothed track
############################################################

############################################################
# LIBRARIES
############################################################

library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)

set.seed(123)

############################################################
# SOURCE MODEL FILES
############################################################

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")
source("simulate_CTCRW_3D.R")

############################################################
# LOAD DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")

############################################################
# DO NOT OVERWRITE HDOP
#
# If your CSV contains a real HDOP column, keep it.
# If it does not, create one.
############################################################

if (!"hdop" %in% names(whale)) {
  whale$hdop <- 0
}

############################################################
# CONVERT TIME
############################################################

whale$time <- ymd_hms(whale$time, quiet = TRUE)

############################################################
# REMOVE ROWS WITH INVALID TIME
############################################################

whale <- whale %>% filter(!is.na(time))

############################################################
# SORT BY TIME
############################################################

whale <- whale %>% arrange(time)

############################################################
# CREATE AUGMENTED DATA
############################################################

aug <- whale %>%
  mutate(
    
    Time_sec =
      as.numeric(
        difftime(
          time,
          min(time),
          units = "secs"
        )
      ),
    
    orig_index =
      seq_len(n())
  )


############################################################
# OBSERVATION MATRIX
############################################################

y <- as.matrix(aug[, c("x", "y", "depth")])

N <- nrow(aug)

############################################################
# INITIAL VALUES
#
# Based on Michelot & Blackwell starting values
############################################################

tau_horiz_hours <- 8
tau_vert_hours <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec <- tau_vert_hours * 3600

beta1_start <- 1 / tau_horiz_sec
beta2_start <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert <- 0.5

sigma1_start <- v_rms_horiz * sqrt(2 * beta1_start)
sigma2_start <- v_rms_vert * sqrt(2 * beta2_start)

params_start <- c(
  beta1 = log(beta1_start),
  beta2 = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

cat("\n")
cat("INITIAL PARAMETERS\n")
cat("==================\n")
print(params_start)

############################################################
# NEGATIVE LOG-LIKELIHOOD
#
# Small measurement error is used to stabilize the model.
############################################################

neg_loglikelihood_noerror <- function(
    params,
    data_aug
) {
  
  ##########################################################
  # Transform parameters
  ##########################################################
  
  beta1 <- exp(params["beta1"])
  beta2 <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  ##########################################################
  # Observations
  ##########################################################
  
  y <- as.matrix(data_aug[, c("x", "y", "depth")])
  N <- nrow(data_aug)
  
  ##########################################################
  # Time
  ##########################################################
  
  time_sec <-
    as.numeric(
      difftime(
        data_aug$time,
        min(data_aug$time),
        units = "secs"
      )
    )
  
  ##########################################################
  # CORRECT TIME INTERVALS
  #
  # First observation is the initial state.
  #
  # delta[1] = 0
  # delta[2] = time[2] - time[1]
  # delta[3] = time[3] - time[2]
  # etc.
  ##########################################################
  
  delta <- c(0,diff(time_sec))
  
  delta[!is.finite(delta)] <- 0
  
  delta <- pmax(delta,0)
  
  ##########################################################
  # Process noise
  ##########################################################
  
  s_horiz <- sigma1^2
  s_vert <- sigma2^2
  
  ##########################################################
  # Parameter vectors
  ##########################################################
  
  beta1_vec <-rep(beta1,N)
  beta2_vec <- rep(beta2,N)
  
  ##########################################################
  # SMALL MEASUREMENT ERROR
  #
  # x = 5 m SD
  # y = 5 m SD
  # depth = 10 m SD
  ##########################################################
  
  Hmat <- matrix(0,N,3)
  Hmat[,1] <- 5^2
  Hmat[,2] <- 5^2
  Hmat[,3] <- 10^2
  
  ##########################################################
  # INITIAL STATE
  ##########################################################
  
  get_first_non_missing <- function(x) {
    idx <- which(!is.na(x))[1]
    if (length(idx) == 0) {
      return(0)
    }
    
    x[idx]
  }
  
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  ##########################################################
  # INITIAL COVARIANCE
  ##########################################################
  
  P <- diag(6) * 1e6
  
  ##########################################################
  # KALMAN FILTER
  ##########################################################
  
  filt <- CTCRW_filter1(
    y = y,
    Hmat = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s_horiz = s_horiz,
    s_vert = s_vert,
    delta = delta,
    a = a,
    P = P
  )
  
  ##########################################################
  # Return negative log-likelihood
  ##########################################################
  
  if (!is.finite(filt$ll)) {
    return(1e100)
  }
  
  -filt$ll
}

############################################################
# FIT MODEL
############################################################

cat("\n")
cat("FITTING CTCRW MODEL...\n")
cat("======================\n")

fit <- optim(
  par = params_start,
  fn = function(p) {neg_loglikelihood_noerror(p, aug)},
  method = "L-BFGS-B",
  control = list(trace = 1, maxit = 1000)
)

############################################################
# CHECK OPTIMIZATION
############################################################

cat("\n")
cat("OPTIMIZATION RESULT\n")
cat("===================\n")

print(fit)

############################################################
# TRANSFORM ESTIMATED PARAMETERS
############################################################

p_hat <- exp(fit$par)

cat("\n")
cat("ESTIMATED PARAMETERS\n")
cat("====================\n")

print(p_hat)

beta1_hat <- p_hat["beta1"]
beta2_hat <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta1_hat
tau2_hat <- 1 / beta2_hat

nu1_hat <- sigma1_hat / sqrt(2 * beta1_hat)
nu2_hat <- sigma2_hat / sqrt(2 * beta2_hat)

cat("\n")
cat("INTERPRETABLE PARAMETERS\n")
cat("========================\n")

cat("Horizontal autocorrelation time tau1 (seconds): ", tau1_hat, "\n")
cat("Horizontal autocorrelation time tau1 (hours):   ", tau1_hat / 3600, "\n")
cat("Vertical autocorrelation time tau2 (seconds):   ", tau2_hat, "\n")
cat("Vertical autocorrelation time tau2 (hours):     ", tau2_hat / 3600, "\n")
cat("Horizontal RMS speed nu1 (m/s):                 ", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):                   ", nu2_hat, "\n")

############################################################
# FILTER + SMOOTHER USING ESTIMATED PARAMETERS
############################################################

############################################################
# Correct delta vector
############################################################

delta <- c(0, diff(aug$Time_sec))

delta[!is.finite(delta)] <- 0

delta <- pmax(delta,0)

############################################################
# Measurement-error matrix
############################################################

Hmat_zero <- matrix(0,N,3)
Hmat_zero[,1] <- 5^2
Hmat_zero[,2] <- 5^2
Hmat_zero[,3] <- 10^2

############################################################
# INITIAL STATE
############################################################

get_first_non_missing <- function(x) {
  idx <- which(!is.na(x))[1]
  if (length(idx) == 0) {
    return(0)
  }
  x[idx]
}


a0 <- c(get_first_non_missing(y[,1]),0,get_first_non_missing(y[,2]),0,
        get_first_non_missing(y[,3]),0
)

############################################################
# INITIAL COVARIANCE
############################################################

P0 <- diag(c(100, 10, 100, 10, 100, 10))

############################################################
# KALMAN FILTER
############################################################

cat("\n")
cat("RUNNING KALMAN FILTER...\n")

filt <- CTCRW_filter1(
  y = y,
  Hmat = Hmat_zero,
  beta1_vec = rep(beta1_hat,N),
  beta2_vec = rep(beta2_hat,N),
  s_horiz = sigma1_hat^2,
  s_vert = sigma2_hat^2,
  delta = delta,
  a = a0,
  P = P0
)

############################################################
# KALMAN SMOOTHER
############################################################

cat("RUNNING KALMAN SMOOTHER...\n")

smooth <- CTCRW_smoother1(
  filter_out = filt,
  beta1_vec = rep(beta1_hat,N),
  beta2_vec = rep(beta2_hat,N),
  s_horiz = sigma1_hat^2,
  s_vert = sigma2_hat^2,
  delta = delta
)

############################################################
# CREATE SMOOTHED TRACK DATA FRAME
############################################################

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time

############################################################
# IMPORTANT:
#
# Do NOT remove rows where the original observation is NA.
#
# The entire point of the smoother is to estimate the latent
# trajectory at times where observations are missing.
############################################################

############################################################
# COMBINE OBSERVATIONS AND SMOOTHED TRACK
############################################################

plot_data <- aug %>%
  mutate(
    smooth_x = smooth_track$x,
    smooth_y = smooth_track$y,
    smooth_depth = smooth_track$depth
  )

############################################################
# FIRST 1000 ROWS
############################################################

n_plot <- min(1000, nrow(plot_data))
subset_plot <- plot_data[1:n_plot, ]

############################################################
# PLOT 1
# DEPTH VS TIME
############################################################

p_depth <- ggplot() +
  
  # Observed depth
  geom_point(data = subset_plot,
    aes(x = time,y = depth),
    color = "blue",
    alpha = 0.5
  ) +
  
  # Smoothed latent depth
  geom_line(data = subset_plot,
    aes(x = time,y = smooth_depth),
    color = "red",
    linewidth = 1
  ) +
  
  labs(title = paste0("Depth vs Time — CTCRW Smoothed Track (First ",n_plot," rows)"),
    subtitle = "Blue = observations; red = smoothed latent trajectory",
    x = "Time",
    y = "Depth"
  ) +
  
  theme_minimal()

print(p_depth)

############################################################
# PLOT 2
# X VS TIME
############################################################

p_x <- ggplot() +
  
  # Observed x
  geom_point(data = subset_plot,
    aes(x = time,y = x),
    color = "blue",
    alpha = 0.5
  ) +
  
  # Smoothed x
  geom_line(data = subset_plot,
    aes(x = time,y = smooth_x),
    color = "red",
    linewidth = 1
  ) +
  
  labs(title = paste0("X vs Time — CTCRW Smoothed Track (First ",n_plot," rows)"),
    subtitle = "Blue = observations; red = smoothed latent trajectory",
    x = "Time",
    y = "X"
  ) +
  
  theme_minimal()

print(p_x)

############################################################
# PLOT 3
# Y VS TIME
############################################################

p_y <- ggplot() +
  
  # Observed y
  geom_point(data = subset_plot,
    aes(x = time,y = y),
    color = "blue",
    alpha = 0.5
  ) +
  
  # Smoothed y
  geom_line(data = subset_plot,
    aes(x = time,y = smooth_y),
    color = "red",
    linewidth = 1
  ) +
  
  labs(title = paste0("Y vs Time — CTCRW Smoothed Track (First ",n_plot," rows)"),
    subtitle = "Blue = observations; red = smoothed latent trajectory",
    x = "Time",
    y = "Y"
  ) +
  
  theme_minimal()

print(p_y)

############################################################
# OPTIONAL: SAVE SMOOTHED TRACK
############################################################

write.csv(
  smooth_track,
  "sperm_whale_CTCRW_smoothed_track.csv",
  row.names = FALSE
)

############################################################
# OPTIONAL SUMMARY
############################################################

cat("\n")
cat("SMOOTHED TRACK CREATED\n")
cat("======================\n")

cat(
  "Number of observations: ",
  N,
  "\n"
)

cat(
  "Number of smoothed states: ",
  nrow(smooth_track),
  "\n"
)

cat(
  "Smoothed track saved to:\n",
  "sperm_whale_CTCRW_smoothed_track.csv\n"
)

cat("\n")
cat("DONE.\n")






