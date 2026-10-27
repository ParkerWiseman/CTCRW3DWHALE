

############################################################
# SpermWhale_2D_CTCRW.R
#
# Fit supervisor's 2D CTCRW model to sperm whale data
############################################################

library(dplyr)
library(lubridate)
library(ggplot2)

############################################################
# SOURCE FILES
############################################################

source("kalman_matrices.R")
source("filter_ctcrw.R")
source("smoother_ctcrw.R")

############################################################
# LOAD DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")

whale$time <- ymd_hms(whale$time, quiet = TRUE)

whale <- whale %>%
  filter(
    !is.na(time),
    !is.na(x),
    !is.na(y)
  ) %>%
  arrange(time)

############################################################
# CREATE TIME AND OBSERVATION MATRICES
############################################################

times <- as.numeric(
  difftime(
    whale$time,
    min(whale$time),
    units = "secs"
  )
)

obs <- as.matrix(
  whale[, c("x","y")]
)

n <- nrow(obs)

############################################################
# MEASUREMENT ERROR
############################################################

H_array <- array(0, c(2,2,n))

for(i in 1:n) {
  H_array[,,i] <- diag(c(5^2, 5^2))
}

############################################################
# NEGATIVE LOG LIKELIHOOD
############################################################

neg_loglikelihood_2d <- function(
    params,
    times,
    obs,
    H_array
) {
  
  tau <- exp(params["log_tau"])
  nu  <- exp(params["log_nu"])
  
  n <- nrow(obs)
  
  tau_vec <- rep(tau, n)
  nu_vec  <- rep(nu,  n)
  
  a0 <- c(
    obs[1,1],
    0,
    obs[1,2],
    0
  )
  
  P0 <- diag(c(100,10,100,10))
  
  filt <- filter_ctcrw(
    times   = times,
    obs     = obs,
    tau     = tau_vec,
    nu      = nu_vec,
    a0      = a0,
    P0      = P0,
    H_array = H_array
  )
  
  ll <- 0
  
  for(i in 1:(n-1)) {
    
    if(any(is.na(obs[i,]))) next
    
    Fi <- filt$F[,,i]
    vi <- matrix(filt$v[i,], ncol = 1)
    
    ll <- ll -
      0.5 * (
        log(det(Fi)) +
          t(vi) %*% solve(Fi) %*% vi +
          2 * log(2*pi)
      )
  }
  
  return(-as.numeric(ll))
}

############################################################
# INITIAL VALUES
############################################################

tau_start <- 8 * 3600
nu_start  <- 1

params_start <- c(
  log_tau = log(tau_start),
  log_nu  = log(nu_start)
)

############################################################
# FIT MODEL
############################################################

cat("\nFITTING 2D CTCRW MODEL\n")
cat("======================\n")

fit <- optim(
  par     = params_start,
  fn      = neg_loglikelihood_2d,
  times   = times,
  obs     = obs,
  H_array = H_array,
  method  = "L-BFGS-B",
  hessian = TRUE,
  control = list(
    trace = 1,
    maxit = 1000
  )
)

############################################################
# PARAMETER ESTIMATES
############################################################

tau_hat <- exp(fit$par["log_tau"])
nu_hat  <- exp(fit$par["log_nu"])

cat("\nESTIMATED PARAMETERS\n")
cat("====================\n")
cat("tau =", tau_hat, "seconds\n")
cat("tau =", tau_hat/3600, "hours\n")
cat("nu  =", nu_hat, "m/s\n")

############################################################
# DELTA-METHOD CIs
############################################################

Sigma_log <- solve(fit$hessian)

se_tau <- tau_hat *
  sqrt(Sigma_log["log_tau","log_tau"])

tau_ci <- c(
  tau_hat - 1.96 * se_tau,
  tau_hat + 1.96 * se_tau
)

se_nu <- nu_hat *
  sqrt(Sigma_log["log_nu","log_nu"])

nu_ci <- c(
  nu_hat - 1.96 * se_nu,
  nu_hat + 1.96 * se_nu
)

cat("\n95% CONFIDENCE INTERVALS\n")
cat("========================\n")

cat(
  "tau: ",
  tau_ci[1],
  " to ",
  tau_ci[2],
  "\n"
)

cat(
  "nu: ",
  nu_ci[1],
  " to ",
  nu_ci[2],
  "\n"
)

############################################################
# FILTER
############################################################

a0 <- c(
  obs[1,1],
  0,
  obs[1,2],
  0
)

P0 <- diag(c(100,10,100,10))

tau_vec <- rep(tau_hat, n)
nu_vec  <- rep(nu_hat, n)

filt <- filter_ctcrw(
  times   = times,
  obs     = obs,
  tau     = tau_vec,
  nu      = nu_vec,
  a0      = a0,
  P0      = P0,
  H_array = H_array
)

############################################################
# SMOOTHER
############################################################

smooth <- smoother_ctcrw(
  times   = times,
  obs     = obs,
  tau     = tau_vec,
  nu      = nu_vec,
  a0      = a0,
  P0      = P0,
  H_array = H_array
)

############################################################
# CREATE SMOOTHED TRACK
############################################################

smooth_track <- data.frame(
  time = whale$time,
  x = smooth$pred[,1],
  vx = smooth$pred[,2],
  y = smooth$pred[,3],
  vy = smooth$pred[,4]
)

n_100 <- min(100, nrow(plot_obs))

plot_obs_100 <- plot_obs[1:n_100, ]
plot_smooth_100 <- plot_smooth[1:n_100, ]


ggplot() +
  geom_point(
    data = plot_obs_100,
    aes(time, x),
    colour = "blue",
    alpha = 0.5
  ) +
  geom_line(
    data = plot_smooth_100,
    aes(time, x),
    colour = "red",
    linewidth = 1
  ) +
  geom_point(
    data = plot_smooth_100,
    aes(time, x),
    colour = "orange",
    size = 1
  ) +
  labs(
    title = paste0(
      "X vs Time (First ",
      n_100,
      " Observations)"
    ),
    subtitle =
      "Blue = observations; Orange = smoothed states; Red = line connecting smoothed states",
    x = "Time",
    y = "X"
  ) +
  theme_minimal()


ggplot() +
  geom_point(
    data = plot_obs_100,
    aes(time, y),
    colour = "blue",
    alpha = 0.5
  ) +
  geom_line(
    data = plot_smooth_100,
    aes(time, y),
    colour = "red",
    linewidth = 1
  ) +
  geom_point(
    data = plot_smooth_100,
    aes(time, y),
    colour = "orange",
    size = 1
  ) +
  labs(
    title = paste0(
      "Y vs Time (First ",
      n_100,
      " Observations)"
    ),
    subtitle =
      "Blue = observations; Orange = smoothed states; Red = line connecting smoothed states",
    x = "Time",
    y = "Y"
  ) +
  theme_minimal()


############################################################
# PLOTS
############################################################

n_plot <- min(1000, nrow(smooth_track))

plot_obs <- whale[1:n_plot, ]
plot_smooth <- smooth_track[1:n_plot, ]

ggplot() +
  geom_point(
    data = plot_obs,
    aes(time, x),
    colour = "blue",
    alpha = 0.5
  ) +
  geom_line(
    data = plot_smooth,
    aes(time, x),
    colour = "red",
    linewidth = 1
  ) +
  geom_point(
    data = plot_smooth,
    aes(time, x),
    colour = "orange",
    size = 1
  ) +
  labs(
    title = "X vs Time",
    subtitle = "Blue = observations; Orange = smoothed states; Red = line connecting smoothed states"
  ) +
  theme_minimal()

ggplot() +
  geom_point(
    data = plot_obs,
    aes(time, y),
    colour = "blue",
    alpha = 0.5
  ) +
  geom_line(
    data = plot_smooth,
    aes(time, y),
    colour = "red",
    linewidth = 1
  ) +
  geom_point(
    data = plot_smooth,
    aes(time, y),
    colour = "orange",
    size = 1
  ) +
  labs(
    title = "Y vs Time",
    subtitle = "Blue = observations; Orange = smoothed states; Red = line connecting smoothed states"
  ) +
  theme_minimal()

############################################################
# SAVE SMOOTHED TRACK
############################################################

write.csv(
  smooth_track,
  "sperm_whale_2D_CTCRW_smoothed_track.csv",
  row.names = FALSE
)

cat("\nDONE\n")











