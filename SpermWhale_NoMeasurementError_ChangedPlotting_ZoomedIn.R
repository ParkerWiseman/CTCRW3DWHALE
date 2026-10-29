


############################################################
# SpermWhale_NoMeasurementError_Good.R
#
# 3D CTCRW sperm whale model (supervisor-style filter + RTS smoother)
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
source("utility.R")

############################################################
# LOAD DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")

############################################################
# DO NOT OVERWRITE HDOP
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
# REPLACE MISSING DEPTHS WITH 0 WHEN POSITION IS OBSERVED
############################################################

whale <- whale %>%
  mutate(
    depth = ifelse(
      !is.na(x) & !is.na(y) & is.na(depth),
      0,
      depth
    )
  )

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
    orig_index = seq_len(n())
  )

############################################################
# OBSERVATION MATRIX
############################################################

y <- as.matrix(aug[, c("x", "y", "depth")])
N <- nrow(aug)

############################################################
# INITIAL VALUES (Michelot & Blackwell style)
############################################################

tau_horiz_hours <- 8
tau_vert_hours  <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours * 3600

beta_horiz_start <- 1 / tau_horiz_sec
beta_vert_start  <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

sigma_horiz_start <- v_rms_horiz * sqrt(2 * beta_horiz_start)
sigma_vert_start  <- v_rms_vert  * sqrt(2 * beta_vert_start)

params_start <- c(
  beta_horiz  = log(beta_horiz_start),
  beta_vert   = log(beta_vert_start),
  sigma_horiz = log(sigma_horiz_start),
  sigma_vert  = log(sigma_vert_start)
)

cat("\nINITIAL PARAMETERS\n==================\n")
print(params_start)



############################################################
# FIT MODEL
############################################################

cat("\nFITTING CTCRW MODEL...\n======================\n")

fit <- optim(
  par     = params_start,
  fn      = function(p) { neg_loglikelihood_noerror(p, aug) },
  method  = "L-BFGS-B",
  hessian = TRUE,
  control = list(trace = 1, maxit = 1000)
)

cat("\nOPTIMIZATION RESULT\n===================\n")
print(fit)


##############################################################################

############################################################
# TRANSFORM ESTIMATED PARAMETERS
############################################################

p_hat <- exp(fit$par)

cat("\nESTIMATED PARAMETERS\n====================\n")
print(p_hat)

beta_horiz_hat  <- p_hat["beta_horiz"]
beta_vert_hat   <- p_hat["beta_vert"]
sigma_horiz_hat <- p_hat["sigma_horiz"]
sigma_vert_hat  <- p_hat["sigma_vert"]

############################################################
# HESSIAN-BASED COVARIANCE MATRIX
############################################################

H <- fit$hessian
Sigma_log <- solve(H)

############################################################
# CONFIDENCE INTERVALS FOR MODEL PARAMETERS
############################################################

se_log <- sqrt(diag(Sigma_log))
se_original <- p_hat * se_log

lower <- p_hat - 1.96 * se_original
upper <- p_hat + 1.96 * se_original

cat("\nPARAMETER CONFIDENCE INTERVALS\n==============================\n")

for(i in seq_along(p_hat)) {
  cat(
    names(p_hat)[i],
    ": ",
    lower[i],
    " to ",
    upper[i],
    "\n"
  )
}

############################################################
# PARAMETER SUMMARY TABLE
############################################################

param_table <- data.frame(
  Parameter = names(p_hat),
  Estimate  = p_hat,
  StdError  = se_original,
  Lower95   = lower,
  Upper95   = upper,
  row.names = NULL
)

cat("\nPARAMETER SUMMARY TABLE\n=======================\n")
print(param_table)

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta_horiz_hat
tau2_hat <- 1 / beta_vert_hat

nu1_hat <- sigma_horiz_hat * sqrt(pi) / (2 * sqrt(beta_horiz_hat))

nu2_hat <- sigma_vert_hat * sqrt(pi) / (2 * sqrt(beta_vert_hat))

############################################################
# DELTA-METHOD STANDARD ERRORS FOR TAU
############################################################

var_log_beta_h <- Sigma_log["beta_horiz", "beta_horiz"]
var_log_beta_v <- Sigma_log["beta_vert",  "beta_vert"]

se_tau1 <- tau1_hat * sqrt(var_log_beta_h)
se_tau2 <- tau2_hat * sqrt(var_log_beta_v)

tau1_lower <- tau1_hat - 1.96 * se_tau1
tau1_upper <- tau1_hat + 1.96 * se_tau1

tau2_lower <- tau2_hat - 1.96 * se_tau2
tau2_upper <- tau2_hat + 1.96 * se_tau2

############################################################
# DELTA-METHOD STANDARD ERRORS FOR NU
############################################################

Sigma_h <- Sigma_log[c("beta_horiz","sigma_horiz"),c("beta_horiz","sigma_horiz")]

grad_nu1 <- c(-nu1_hat / 2, nu1_hat)

se_nu1 <- sqrt(as.numeric(t(grad_nu1) %*% Sigma_h %*% grad_nu1))

nu1_lower <- nu1_hat - 1.96 * se_nu1
nu1_upper <- nu1_hat + 1.96 * se_nu1

Sigma_v <- Sigma_log[c("beta_vert","sigma_vert"),c("beta_vert","sigma_vert")]

grad_nu2 <- c(-nu2_hat / 2, nu2_hat)

se_nu2 <- sqrt(as.numeric(t(grad_nu2) %*% Sigma_v %*% grad_nu2))

nu2_lower <- nu2_hat - 1.96 * se_nu2
nu2_upper <- nu2_hat + 1.96 * se_nu2

############################################################
# PRINT INTERPRETABLE PARAMETERS
############################################################

cat("\nINTERPRETABLE PARAMETERS\n========================\n")

cat("Horizontal autocorrelation time tau1 (seconds): ",tau1_hat,"\n")

cat("Horizontal autocorrelation time tau1 (hours): ",tau1_hat / 3600,"\n")

cat("Vertical autocorrelation time tau2 (seconds): ",tau2_hat,"\n")

cat("Vertical autocorrelation time tau2 (hours): ",tau2_hat / 3600,"\n")

cat("Horizontal RMS speed nu1 (m/s): ",nu1_hat,"\n")

cat("Vertical RMS speed nu2 (m/s): ",nu2_hat,"\n")

############################################################
# DELTA-METHOD CONFIDENCE INTERVALS
############################################################

interp_table <- data.frame(
  Parameter = c(
    "tau1_sec",
    "tau2_sec",
    "nu1",
    "nu2"
  ),
  Estimate = c(
    tau1_hat,
    tau2_hat,
    nu1_hat,
    nu2_hat
  ),
  StdError = c(
    se_tau1,
    se_tau2,
    se_nu1,
    se_nu2
  ),
  Lower95 = c(
    tau1_lower,
    tau2_lower,
    nu1_lower,
    nu2_lower
  ),
  Upper95 = c(
    tau1_upper,
    tau2_upper,
    nu1_upper,
    nu2_upper
  )
)

cat(
  "\nDELTA-METHOD CONFIDENCE INTERVALS\n",
  "=================================\n",
  sep = ""
)

print(interp_table)






#############################################################################


############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta_horiz_hat
tau2_hat <- 1 / beta_vert_hat

nu1_hat <- sigma_horiz_hat * sqrt(pi) / (2 * sqrt(beta_horiz_hat))
nu2_hat <- sigma_vert_hat * sqrt(pi)  / (2 * sqrt(beta_vert_hat))

cat("\nINTERPRETABLE PARAMETERS\n========================\n")
cat("Horizontal autocorrelation time tau1 (seconds): ", tau1_hat, "\n")
cat("Horizontal autocorrelation time tau1 (hours):   ", tau1_hat / 3600, "\n")
cat("Vertical autocorrelation time tau2 (seconds):   ", tau2_hat, "\n")
cat("Vertical autocorrelation time tau2 (hours):     ", tau2_hat / 3600, "\n")
cat("Horizontal RMS speed nu1 (m/s):                 ", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):                   ", nu2_hat, "\n")


############################################################
# FILTER + SMOOTHER USING ESTIMATED PARAMETERS
############################################################

delta <- c(diff(aug$Time_sec),0)


Hmat_zero <- matrix(0, N, 3)
Hmat_zero[,1] <- 100^2
Hmat_zero[,2] <- 100^2
Hmat_zero[,3] <- 10^2


a1 <- c(
  get_first_non_missing(y[,1]), 0,
  get_first_non_missing(y[,2]), 0,
  get_first_non_missing(y[,3]), 0
)


P1 <- diag(c(100, 10, 100, 10, 100, 10))

cat("\nRUNNING KALMAN FILTER...\n")

filt <- CTCRW_filter(
  y              = y,
  Hmat           = Hmat_zero,
  beta_horiz_vec = rep(beta_horiz_hat,  N),
  beta_vert_vec  = rep(beta_vert_hat,   N),
  var_horiz_vec  = rep(sigma_horiz_hat^2, N),
  var_vert_vec   = rep(sigma_vert_hat^2,  N),
  delta          = delta,
  a              = a1,
  P              = P1
)

cat("RUNNING KALMAN SMOOTHER...\n")

smooth <- CTCRW_smoother(filt)

############################################################
# CREATE SMOOTHED TRACK DATA FRAME
############################################################

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time

############################################################
# COMBINE OBSERVATIONS AND SMOOTHED TRACK
############################################################

plot_data <- aug %>%
  mutate(
    smooth_x     = smooth_track$x,
    smooth_y     = smooth_track$y,
    smooth_depth = smooth_track$depth
  )


############################################################
# FIRST 100 HORIZONTAL OBSERVATIONS
############################################################

xy_plot <- plot_data %>%
  filter(
    !is.na(x),
    !is.na(y)
  )

n_xy <- min(100, nrow(xy_plot))

subset_xy <- xy_plot[1:n_xy, ]



n_plot <- min(1000, nrow(plot_data))
subset_plot <- plot_data[1:n_plot, ]


############################################################
# FIRST 200 OBSERVATIONS
############################################################

n_200 <- min(200, nrow(plot_data))
subset_200 <- plot_data[1:n_200, ]



############################################################
# PLOT 1.0: DEPTH VS TIME (FIRST 200 OBSERVATIONS)
############################################################

p_depth_200 <- ggplot() +
  
  # green line connecting observations
  geom_line(
    data = subset_200,
    aes(x = time, y = depth),
    color = "green",
    linewidth = 1
  ) +
  
  # red line connecting smoothed states
  geom_line(
    data = subset_200,
    aes(x = time, y = smooth_depth),
    color = "red",
    linewidth = 1
  ) +
  
  # observed depths
  geom_point(
    data = subset_200,
    aes(x = time, y = depth),
    color = "blue",
    alpha = 0.5
  ) +
  
  # orange smoothed states
  geom_point(
    data = subset_200,
    aes(x = time, y = smooth_depth),
    color = "orange",
    size = 1
  ) +
  
  labs(
    title = paste0(
      "Depth vs Time — CTCRW Smoothed Track (First ",
      n_200,
      " rows)"
    ),
    subtitle =
      "Blue = observations; Green = line segments connecting observations; Orange = smoothed states; Red = line segments connecting smoothed states",
    x = "Time",
    y = "Depth"
  ) +
  
  theme_minimal()

print(p_depth_200)






############################################################
# PLOT 1: DEPTH VS TIME
############################################################

p_depth <- ggplot() +
  
  # green line connecting observations
  geom_line(
    data = subset_plot,
    aes(x = time, y = depth),
    color = "green",
    linewidth = 1
  ) +
  
  # red line connecting smoothed states
  geom_line(
    data = subset_plot,
    aes(x = time, y = smooth_depth),
    color = "red",
    linewidth = 1
  ) +
  
  # observed depths
  geom_point(
    data = subset_plot,
    aes(x = time, y = depth),
    color = "blue",
    alpha = 0.5
  ) +
  
  # orange smoothed states
  geom_point(
    data = subset_plot,
    aes(x = time, y = smooth_depth),
    color = "orange",
    size = 1
  ) +
  
  labs(
    title = paste0(
      "Depth vs Time — CTCRW Smoothed Track (First ",
      n_plot,
      " rows)"
    ),
    subtitle =
      "Blue = observations; Green = line segments connecting observations; Orange = smoothed states; Red = line segments connecting smoothed states",
    x = "Time",
    y = "Depth"
  ) +
  
  theme_minimal()

print(p_depth)



############################################################
# PLOT 2.0: X VS TIME (FIRST 100 X OBSERVATIONS)
############################################################

p_x <- ggplot() +
  
  # green line connecting observations
  geom_line(
    data = subset_xy,
    aes(x = time, y = x),
    color = "green",
    linewidth = 1
  ) +
  
  # red line connecting smoothed states
  geom_line(
    data = subset_xy,
    aes(x = time, y = smooth_x),
    color = "red",
    linewidth = 1
  ) +
  
  geom_point(
    data = subset_xy,
    aes(x = time, y = x),
    color = "blue",
    alpha = 0.5
  ) +
  
  geom_point(
    data = subset_xy,
    aes(x = time, y = smooth_x),
    color = "orange",
    size = 1
  ) +
  
  labs(
    title = paste0(
      "X vs Time — First ",
      n_xy,
      " X Observations"
    ),
    subtitle =
      "Blue = observations; Green = line segments connecting observations; Orange = smoothed states; Red = line segments connecting smoothed states",
    x = "Time",
    y = "X"
  ) +
  
  theme_minimal()

print(p_x)







############################################################
# PLOT 2: X VS TIME
############################################################

p_x <- ggplot() +
  
  # green line connecting observations
  geom_line(
    data = subset_plot,
    aes(x = time, y = x),
    color = "green",
    linewidth = 1
  ) +
  
  # red line connecting smoothed states
  geom_line(
    data = subset_plot,
    aes(x = time, y = smooth_x),
    color = "red",
    linewidth = 1
  ) +
  
  # observed positions
  geom_point(
    data = subset_plot,
    aes(x = time, y = x),
    color = "blue",
    alpha = 0.5
  ) +
  
  # orange smoothed states
  geom_point(
    data = subset_plot,
    aes(x = time, y = smooth_x),
    color = "orange",
    size = 1
  ) +
  
  labs(
    title = paste0(
      "X vs Time — CTCRW Smoothed Track (First ",
      n_plot,
      " rows)"
    ),
    subtitle =
      "Blue = observations; Green = line segments connecting observations; Orange = smoothed states; Red = line segments connecting smoothed states",
    x = "Time",
    y = "X"
  ) +
  
  theme_minimal()

print(p_x)


############################################################
# PLOT 3.0: Y VS TIME (FIRST 100 Y OBSERVATIONS)
############################################################

p_y <- ggplot() +
  
  # green line connecting observations
  geom_line(
    data = subset_xy,
    aes(x = time, y = y),
    color = "green",
    linewidth = 1
  ) +
  
  # red line connecting smoothed states
  geom_line(
    data = subset_xy,
    aes(x = time, y = smooth_y),
    color = "red",
    linewidth = 1
  ) +
  
  geom_point(
    data = subset_xy,
    aes(x = time, y = y),
    color = "blue",
    alpha = 0.5
  ) +
  
  geom_point(
    data = subset_xy,
    aes(x = time, y = smooth_y),
    color = "orange",
    size = 1
  ) +
  
  labs(
    title = paste0(
      "Y vs Time — First ",
      n_xy,
      " Y Observations"
    ),
    subtitle =
      "Blue = observations; Green = line segments connecting observations; Orange = smoothed states; Red = line segments connecting smoothed states",
    x = "Time",
    y = "Y"
  ) +
  
  theme_minimal()

print(p_y)


############################################################
# PLOT 3: Y VS TIME
############################################################

p_y <- ggplot() +
  
  # green line connecting observations
  geom_line(
    data = subset_plot,
    aes(x = time, y = y),
    color = "green",
    linewidth = 1
  ) +
  
  # red line connecting smoothed states
  geom_line(
    data = subset_plot,
    aes(x = time, y = smooth_y),
    color = "red",
    linewidth = 1
  ) +
  
  # observed positions
  geom_point(
    data = subset_plot,
    aes(x = time, y = y),
    color = "blue",
    alpha = 0.5
  ) +
  
  # orange smoothed states
  geom_point(
    data = subset_plot,
    aes(x = time, y = smooth_y),
    color = "orange",
    size = 1
  ) +
  
  labs(
    title = paste0(
      "Y vs Time — CTCRW Smoothed Track (First ",
      n_plot,
      " rows)"
    ),
    subtitle =
      "Blue = observations; Green = line segments connecting observations Orange = smoothed states; Red = line segments connecting smoothed states",
    x = "Time",
    y = "Y"
  ) +
  
  theme_minimal()

print(p_y)

############################################################
# OPTIONAL: SAVE SMOOTHED TRACK
############################################################

write.csv(smooth_track, "sperm_whale_CTCRW_smoothed_track.csv", row.names = FALSE)

cat("\nSMOOTHED TRACK CREATED\n======================\n")
cat("Number of observations: ", N, "\n")
cat("Number of smoothed states: ", nrow(smooth_track), "\n")
cat("Smoothed track saved to: sperm_whale_CTCRW_smoothed_track.csv\n")
cat("\nDONE.\n")













