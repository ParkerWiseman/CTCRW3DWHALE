


############################################################
# SpermWhale_AnalysisCode_WithUnobservedTimePoints_gamma_TimeOfDay.R
#
# Uses time of day (more biologically meaningful)
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
source("CTCRW_filter_gamma.R")
source("CTCRW_smoother_gamma.R")
source("neg_loglikelihood_gamma.R")
source("simulate_CTCRW_3D.R")
source("utility.R")
source("gamma_spline.R")

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

prediction_interval_sec <- 60

regular_grid <- seq(
  min(whale$time),
  max(whale$time),
  by = prediction_interval_sec
)

all_times <- sort(
  unique(
    c(
      whale$time,
      regular_grid
    )
  )
)

aug <- data.frame(
  time = all_times
) %>%
  left_join(
    whale,
    by = "time"
  ) %>%
  arrange(time) %>%
  mutate(
    Time_sec =
      as.numeric(
        difftime(
          time,
          min(time),
          units = "secs"
        )
      )
  )




aug <- aug %>%
  mutate(
    observed =
      !is.na(x) |
      !is.na(y) |
      !is.na(depth)
  )

# rescaling time before constructing the spline, so the spline is
# working in hours rather than seconds.
# This gives much better numerical stability for optimization.
#aug$covariate <- aug$Time_sec / 3600
aug$covariate <- hour(aug$time) + minute(aug$time)/60 + second(aug$time)/3600

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
  sigma_vert  = log(sigma_vert_start),
  
  gamma_intercept = 0,
  gamma_b1 = 0,
  gamma_b2 = 0,
  gamma_b3 = 0,
  gamma_b4 = 0,
  gamma_b5 = 0,
  gamma_b6 = 0
)

cat("\nINITIAL PARAMETERS\n==================\n")
print(params_start)



############################################################
# FIT MODEL
############################################################

cat("\nFITTING CTCRW MODEL...\n======================\n")

fit <- optim(
  par     = params_start,
  fn      = function(p) { neg_loglikelihood_gamma(p, aug) },
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

beta_horiz_hat  <- exp(fit$par["beta_horiz"])
beta_vert_hat   <- exp(fit$par["beta_vert"])

sigma_horiz_hat <- exp(fit$par["sigma_horiz"])
sigma_vert_hat  <- exp(fit$par["sigma_vert"])

gamma_hat <- fit$par[
  c(
    "gamma_intercept",
    "gamma_b1",
    "gamma_b2",
    "gamma_b3",
    "gamma_b4",
    "gamma_b5",
    "gamma_b6"
  )
]


############################################################

gamma_coef_hat <- c(
  fit$par["gamma_b1"],
  fit$par["gamma_b2"],
  fit$par["gamma_b3"],
  fit$par["gamma_b4"],
  fit$par["gamma_b5"],
  fit$par["gamma_b6"]
)

gamma_vec_hat <- make_gamma(
  covariate = aug$covariate,
  gamma_intercept = fit$par["gamma_intercept"],
  gamma_coef = gamma_coef_hat
)

summary(gamma_vec_hat)


cat(
  "\nRange of gamma_t:\n",
  min(gamma_vec_hat),
  "to",
  max(gamma_vec_hat),
  "\n"
)


cat("\nESTIMATED PARAMETERS\n====================\n")

cat("beta_horiz  =", beta_horiz_hat, "\n")
cat("beta_vert   =", beta_vert_hat, "\n")
cat("sigma_horiz =", sigma_horiz_hat, "\n")
cat("sigma_vert  =", sigma_vert_hat, "\n")

cat("\nGamma spline coefficients\n")
print(gamma_hat)

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta_horiz_hat
tau2_hat <- 1 / beta_vert_hat

nu1_hat <- sigma_horiz_hat * sqrt(pi) / (2 * sqrt(beta_horiz_hat))
nu2_hat <- sigma_vert_hat * sqrt(pi) / (2 * sqrt(beta_vert_hat))



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







#############################################################################


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

filt <- CTCRW_filter_gamma(
  y              = y,
  Hmat           = Hmat_zero,
  beta_horiz_vec = rep(beta_horiz_hat,  N),
  beta_vert_vec  = rep(beta_vert_hat,   N),
  var_horiz_vec  = rep(sigma_horiz_hat^2, N),
  var_vert_vec   = rep(sigma_vert_hat^2,  N),
  delta          = delta,
  a              = a1,
  P              = P1,
  gamma_vec      = gamma_vec_hat
)

cat("RUNNING KALMAN SMOOTHER...\n")

smooth <- CTCRW_smoother_gamma(filt)

############################################################
# CREATE SMOOTHED TRACK DATA FRAME
############################################################

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time


smooth_track$gamma <- gamma_vec_hat

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
# PLOT 4: X VS Y
############################################################

xy_data <- plot_data %>%
  filter(
    !is.na(x),
    !is.na(y)
  )

n_xy_plot <- min(1000, nrow(xy_data))

xy_subset <- xy_data[1:n_xy_plot, ]

p_xy <- ggplot() +
  
  # observed track
  geom_path(
    data = xy_subset,
    aes(x = x, y = y),
    color = "green",
    linewidth = 1
  ) +
  
  # smoothed track
  geom_path(
    data = xy_subset,
    aes(x = smooth_x, y = smooth_y),
    color = "red",
    linewidth = 1
  ) +
  
  # observed positions
  geom_point(
    data = xy_subset,
    aes(x = x, y = y),
    color = "blue",
    alpha = 0.5
  ) +
  
  # smoothed positions
  geom_point(
    data = xy_subset,
    aes(x = smooth_x, y = smooth_y),
    color = "orange",
    size = 1
  ) +
  
  labs(
    title = paste0(
      "X vs Y — CTCRW Smoothed Track (First ",
      n_xy_plot,
      " horizontal observations)"
    ),
    subtitle =
      "Blue = observations; Green = observed track; Orange = smoothed states; Red = smoothed track",
    x = "X Position",
    y = "Y Position"
  ) +
  
  coord_equal() +
  theme_minimal()

print(p_xy)


gamma_df <- data.frame(
  time = aug$time,
  gamma = gamma_vec_hat
)

p_gamma <- ggplot(
  gamma_df,
  aes(x = time, y = gamma)
) +
  geom_line() +
  theme_minimal() +
  labs(
    title = expression(hat(gamma)[t]),
    x = "Time",
    y = expression(gamma[t])
  )

print(p_gamma)

############################################################
# OPTIONAL: SAVE SMOOTHED TRACK
############################################################

write.csv(smooth_track, "sperm_whale_CTCRW_smoothed_track.csv", row.names = FALSE)

cat("\nSMOOTHED TRACK CREATED\n======================\n")
cat("Number of observations: ", N, "\n")
cat("Number of smoothed states: ", nrow(smooth_track), "\n")
cat("Smoothed track saved to: sperm_whale_CTCRW_smoothed_track.csv\n")
cat("\nDONE.\n")


ggplot(
  smooth_track,
  aes(x = time)
) +
  geom_line(
    aes(y = depth),
    colour = "blue"
  ) +
  geom_line(
    aes(y = 500 * gamma),
    colour = "red"
  ) +
  theme_minimal()




gamma_df <- data.frame(
  hour_of_day = aug$covariate,
  gamma = gamma_vec_hat
)

ggplot(
  gamma_df,
  aes(
    x = hour_of_day,
    y = gamma
  )
) +
  geom_point(alpha = 0.05) +
  geom_smooth() +
  theme_minimal() +
  labs(
    x = "Hour of day",
    y = expression(gamma[t]),
    title = expression(hat(gamma)[t]~"vs Hour of Day")
  )



summary(smooth_track$vdepth)
summary(gamma_vec_hat)









