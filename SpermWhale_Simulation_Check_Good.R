


############################################################
# SIMULATED SPERM WHALE CTCRW DATA
#
# Purpose:
#
#   1. Generate 5000 observations from known CTCRW
#      parameters.
#
#   2. Create an observation pattern resembling the
#      sperm whale dataset:
#
#        - majority = depth only
#        - minority = horizontal x/y only
#
#   3. Save the resulting dataset.
#
#   4. The resulting CSV can then be supplied directly
#      to the custom CTCRW fitting code.
#
############################################################


############################################################
# 1. LIBRARIES
############################################################

library(MASS)
library(lubridate)
library(dplyr)


############################################################
# 2. SOURCE MODEL FUNCTIONS
############################################################

source("matrices.R")

source("simulate_CTCRW_3D.R")


############################################################
# 3. REPRODUCIBILITY
############################################################

set.seed(12345)


############################################################
# 4. SIMULATION SETTINGS
############################################################

N <- 5000


############################################################
# TIME STEP
############################################################
#
# 60 seconds between latent states.
#
# simulate_CTCRW_3D expects dt in DAYS because the
# function later does:
#
#   time_vec * 86400
#
############################################################

dt_seconds <- 60

dt <- dt_seconds


############################################################
# 5. TRUE CTCRW PARAMETERS
############################################################
#
# These are the values we KNOW are correct because we
# specify them ourselves.
#
############################################################

beta1_true <- 0.0005

beta2_true <- 0.004

sigma1_true <- 0.040

sigma2_true <- 0.080


############################################################
# 6. TRUE INTERPRETABLE PARAMETERS
############################################################

tau1_true <- 1 / beta1_true

tau2_true <- 1 / beta2_true

nu1_true <-
  sigma1_true /
  sqrt(
    2 * beta1_true
  )

nu2_true <-
  sigma2_true /
  sqrt(
    2 * beta2_true
  )


############################################################
# 7. PRINT TRUE PARAMETERS
############################################################

cat("\n")
cat("====================================================\n")
cat("TRUE SIMULATION PARAMETERS\n")
cat("====================================================\n")

cat(
  "beta1:",
  beta1_true,
  "\n"
)

cat(
  "beta2:",
  beta2_true,
  "\n"
)

cat(
  "sigma1:",
  sigma1_true,
  "\n"
)

cat(
  "sigma2:",
  sigma2_true,
  "\n"
)

cat("\n")

cat(
  "Horizontal tau (seconds):",
  tau1_true,
  "\n"
)

cat(
  "Horizontal tau (hours):",
  tau1_true / 3600,
  "\n"
)

cat(
  "Vertical tau (seconds):",
  tau2_true,
  "\n"
)

cat(
  "Vertical tau (hours):",
  tau2_true / 3600,
  "\n"
)

cat("\n")

cat(
  "Horizontal RMS speed (m/s):",
  nu1_true,
  "\n"
)

cat(
  "Vertical RMS speed (m/s):",
  nu2_true,
  "\n"
)


############################################################
# 8. SIMULATE LATENT CTCRW TRAJECTORY
############################################################

cat("\n")
cat("====================================================\n")
cat("SIMULATING LATENT CTCRW TRAJECTORY\n")
cat("====================================================\n")

latent <- simulate_CTCRW_3D(
  
  N = N,
  
  dt = dt,
  
  beta1_true = beta1_true,
  
  beta2_true = beta2_true,
  
  sigma1_true = sigma1_true,
  
  sigma2_true = sigma2_true,
  
  x0 = 0,
  
  y0 = 0,
  
  depth0 = -50
  
)


############################################################
# 9. CHECK SIMULATION
############################################################

cat(
  "Number of simulated states:",
  nrow(latent),
  "\n"
)

cat(
  "Time step (seconds):",
  dt_seconds,
  "\n"
)

cat(
  "Total simulated duration (hours):",
  (N - 1) * dt_seconds / 3600,
  "\n"
)


############################################################
# 10. CREATE OBSERVATION PATTERN
############################################################
#
# We want approximately the same pattern as the real
# sperm whale dataset:
#
#   total = 4468
#   horizontal positions = 731
#   depth-only = 3737
#
# For 5000 simulated observations, we'll use:
#
#   horizontal = 16.4%
#   depth-only  = 83.6%
#
# This is approximately the same proportion as the real
# dataset.
#
############################################################

horizontal_fraction <- 731 / 4468

n_horizontal <- round(
  N * horizontal_fraction
)

n_depth <- N - n_horizontal


############################################################
# 11. RANDOMLY SELECT HORIZONTAL OBSERVATIONS
############################################################
#
# These observations will contain:
#
#   x
#   y
#
# but NOT depth.
#
############################################################

horizontal_index <- sample(
  1:N,
  size = n_horizontal,
  replace = FALSE
)


############################################################
# 12. INITIALIZE OBSERVATION DATASET
############################################################

sim_data <- data.frame(
  
  time = latent$time,
  
  x = NA_real_,
  
  y = NA_real_,
  
  depth = NA_real_
  
)


############################################################
# 13. INSERT HORIZONTAL OBSERVATIONS
############################################################

sim_data$x[
  horizontal_index
] <- latent$x[
  horizontal_index
]

sim_data$y[
  horizontal_index
] <- latent$y[
  horizontal_index
]


############################################################
# 14. IDENTIFY DEPTH OBSERVATIONS
############################################################

depth_index <- setdiff(
  1:N,
  horizontal_index
)


############################################################
# 15. INSERT DEPTH OBSERVATIONS
############################################################

sim_data$depth[
  depth_index
] <- latent$depth[
  depth_index
]


############################################################
# 16. VERIFY OBSERVATION STRUCTURE
############################################################

has_xy <-
  !is.na(sim_data$x) &
  !is.na(sim_data$y)

has_depth <-
  !is.na(sim_data$depth)


############################################################
# 17. PRINT OBSERVATION COUNTS
############################################################

cat("\n")
cat("====================================================\n")
cat("SIMULATED OBSERVATION STRUCTURE\n")
cat("====================================================\n")

cat(
  "Total observations:",
  nrow(sim_data),
  "\n"
)

cat(
  "Horizontal x/y observations:",
  sum(has_xy),
  "\n"
)

cat(
  "Depth observations:",
  sum(has_depth),
  "\n"
)

cat(
  "Depth-only observations:",
  sum(
    has_depth &
      !has_xy
  ),
  "\n"
)

cat(
  "Horizontal-only observations:",
  sum(
    has_xy &
      !has_depth
  ),
  "\n"
)

cat(
  "Observations containing both:",
  sum(
    has_xy &
      has_depth
  ),
  "\n"
)

cat(
  "Completely missing observations:",
  sum(
    !has_xy &
      !has_depth
  ),
  "\n"
)


############################################################
# 18. CHECK THAT THE STRUCTURE IS CORRECT
############################################################

if (
  sum(
    has_xy &
    has_depth
  ) > 0
) {
  
  stop(
    "ERROR: Some observations contain both horizontal and depth measurements."
  )
  
}


if (
  sum(
    !has_xy &
    !has_depth
  ) > 0
) {
  
  stop(
    "ERROR: Some observations contain neither horizontal nor depth measurements."
  )
  
}


############################################################
# 19. ADD TRUE LATENT VALUES
############################################################
#
# These columns are NOT used by the fitting code.
#
# They are included so we can later compare:
#
#   true latent trajectory
#   observed trajectory
#   estimated smoothed trajectory
#
############################################################

sim_data$true_x <- latent$x

sim_data$true_y <- latent$y

sim_data$true_depth <- latent$depth


############################################################
# 20. SAVE SIMULATED DATA
############################################################
#
# This is the file that will be used by your CTCRW model.
#
############################################################

write.csv(
  
  sim_data,
  
  "sperm_whale_simulated_5000.csv",
  
  row.names = FALSE
  
)


############################################################
# 21. SAVE TRUE PARAMETERS
############################################################

true_parameters <- data.frame(
  
  parameter = c(
    "beta1",
    "beta2",
    "sigma1",
    "sigma2",
    "tau1_seconds",
    "tau2_seconds",
    "tau1_hours",
    "tau2_hours",
    "RMS_horizontal",
    "RMS_vertical"
  ),
  
  true_value = c(
    beta1_true,
    beta2_true,
    sigma1_true,
    sigma2_true,
    tau1_true,
    tau2_true,
    tau1_true / 3600,
    tau2_true / 3600,
    nu1_true,
    nu2_true
  )
  
)


write.csv(
  
  true_parameters,
  
  "sperm_whale_simulated_true_parameters.csv",
  
  row.names = FALSE
  
)


############################################################
# 22. FINAL SUMMARY
############################################################

cat("\n")
cat("====================================================\n")
cat("SIMULATION COMPLETE\n")
cat("====================================================\n")

cat(
  "Simulated dataset:",
  "sperm_whale_simulated_5000.csv",
  "\n"
)

cat(
  "True parameters:",
  "sperm_whale_simulated_true_parameters.csv",
  "\n"
)

cat("\n")

cat(
  "Total observations:",
  N,
  "\n"
)

cat(
  "Horizontal-only observations:",
  n_horizontal,
  "\n"
)

cat(
  "Depth-only observations:",
  n_depth,
  "\n"
)

cat("\n")

cat(
  "TRUE beta1:",
  beta1_true,
  "\n"
)

cat(
  "TRUE beta2:",
  beta2_true,
  "\n"
)

cat(
  "TRUE sigma1:",
  sigma1_true,
  "\n"
)

cat(
  "TRUE sigma2:",
  sigma2_true,
  "\n"
)

cat("\n")
cat("DONE.\n")







############################
############################
############################







############################################################
# SpermWhale_Simulated_CTCRW_Validation.R
#
# 3D CTCRW simulation-recovery test
#
# PURPOSE:
#
#   1. Load simulated CTCRW data
#   2. Estimate:
#        beta1
#        beta2
#        sigma1
#        sigma2
#   3. Perform Kalman filtering
#   4. Perform Kalman smoothing
#   5. Plot observations + smoothed track
#   6. Compare estimated parameters with the TRUE
#      parameters used to generate the simulation
#
# IMPORTANT:
#
# The simulated data were generated with NO measurement error.
# Therefore Hmat = 0 is used throughout this validation.
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
# LOAD SIMULATED DATA
############################################################

cat("\n")
cat("LOADING SIMULATED DATA\n")
cat("======================\n")

whale <- read.csv(
  "sperm_whale_simulated_5000.csv"
)

cat(
  "Dataset loaded:",
  nrow(whale),
  "rows\n"
)


############################################################
# LOAD TRUE PARAMETERS
#
# These were saved by the simulation script.
############################################################

true_parameters_file <-
  "sperm_whale_simulated_true_parameters.csv"

if (file.exists(true_parameters_file)) {
  
  true_parameters <- read.csv(
    true_parameters_file
  )
  
  cat("\n")
  cat("TRUE PARAMETERS FILE FOUND\n")
  cat("==========================\n")
  
  print(true_parameters)
  
} else {
  
  cat("\n")
  cat("WARNING: TRUE PARAMETERS FILE NOT FOUND.\n")
  cat("The model will still run, but the final\n")
  cat("true-vs-estimated comparison will not be available.\n")
  
}


############################################################
# CONVERT TIME
############################################################

whale$time <- ymd_hms(
  whale$time,
  quiet = TRUE
)


############################################################
# REMOVE ROWS WITH INVALID TIME
############################################################

whale <- whale %>%
  filter(!is.na(time))


############################################################
# SORT BY TIME
############################################################

whale <- whale %>%
  arrange(time)


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
# CHECK DATA STRUCTURE
############################################################

cat("\n")
cat("SIMULATED DATA STRUCTURE\n")
cat("========================\n")

cat(
  "Total observations:",
  nrow(aug),
  "\n"
)

cat(
  "Horizontal x/y observations:",
  sum(
    !is.na(aug$x) &
      !is.na(aug$y)
  ),
  "\n"
)

cat(
  "Depth observations:",
  sum(!is.na(aug$depth)),
  "\n"
)

cat(
  "Depth-only observations:",
  sum(
    !is.na(aug$depth) &
      is.na(aug$x) &
      is.na(aug$y)
  ),
  "\n"
)

cat(
  "Horizontal-only observations:",
  sum(
    !is.na(aug$x) &
      !is.na(aug$y) &
      is.na(aug$depth)
  ),
  "\n"
)

cat(
  "Observations containing both:",
  sum(
    !is.na(aug$x) &
      !is.na(aug$y) &
      !is.na(aug$depth)
  ),
  "\n"
)

cat(
  "Completely missing observations:",
  sum(
    is.na(aug$x) &
      is.na(aug$y) &
      is.na(aug$depth)
  ),
  "\n"
)


############################################################
# OBSERVATION MATRIX
############################################################

y <- as.matrix(
  aug[, c("x", "y", "depth")]
)

N <- nrow(aug)


############################################################
# INITIAL VALUES
#
# These are ONLY starting values for optimization.
#
# They are NOT the true parameters.
############################################################

tau_horiz_hours <- 8
tau_vert_hours <- 3

tau_horiz_sec <-
  tau_horiz_hours * 3600

tau_vert_sec <-
  tau_vert_hours * 3600

beta1_start <-
  1 / tau_horiz_sec

beta2_start <-
  1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert <- 0.5

sigma1_start <-
  v_rms_horiz *
  sqrt(
    2 * beta1_start
  )

sigma2_start <-
  v_rms_vert *
  sqrt(
    2 * beta2_start
  )

params_start <- c(
  beta1 = log(beta1_start),
  beta2 = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)


cat("\n")
cat("INITIAL OPTIMIZATION PARAMETERS\n")
cat("===============================\n")

cat(
  "beta1 start:",
  beta1_start,
  "\n"
)

cat(
  "beta2 start:",
  beta2_start,
  "\n"
)

cat(
  "sigma1 start:",
  sigma1_start,
  "\n"
)

cat(
  "sigma2 start:",
  sigma2_start,
  "\n"
)


############################################################
# NEGATIVE LOG-LIKELIHOOD
#
# IMPORTANT:
#
# The simulated dataset has NO measurement error.
#
# Therefore:
#
#       Hmat = 0
#
# for x, y, and depth.
############################################################

neg_loglikelihood_noerror <- function(
    params,
    data_aug
) {
  
  
  ##########################################################
  # TRANSFORM PARAMETERS
  ##########################################################
  
  beta1 <- exp(
    params["beta1"]
  )
  
  beta2 <- exp(
    params["beta2"]
  )
  
  sigma1 <- exp(
    params["sigma1"]
  )
  
  sigma2 <- exp(
    params["sigma2"]
  )
  
  
  ##########################################################
  # OBSERVATIONS
  ##########################################################
  
  y <- as.matrix(
    data_aug[, c(
      "x",
      "y",
      "depth"
    )]
  )
  
  N <- nrow(data_aug)
  
  
  ##########################################################
  # TIME
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
  # TIME INTERVALS
  #
  # delta[1] = 0
  #
  # delta[i] =
  # time[i] - time[i-1]
  ##########################################################
  
  delta <- c(
    0,
    diff(time_sec)
  )
  
  delta[!is.finite(delta)] <- 0
  
  delta <- pmax(
    delta,
    0
  )
  
  
  ##########################################################
  # PROCESS NOISE
  ##########################################################
  
  s_horiz <- sigma1^2
  
  s_vert <- sigma2^2
  
  
  ##########################################################
  # PARAMETER VECTORS
  ##########################################################
  
  beta1_vec <- rep(
    beta1,
    N
  )
  
  beta2_vec <- rep(
    beta2,
    N
  )
  
  
  ##########################################################
  # ZERO MEASUREMENT ERROR
  #
  # The simulator generates exact latent observations.
  #
  # Therefore:
  #
  # x measurement variance = 0
  # y measurement variance = 0
  # depth measurement variance = 0
  ##########################################################
  
  Hmat <- matrix(
    0,
    N,
    3
  )
  
  
  ##########################################################
  # INITIAL STATE
  ##########################################################
  
  get_first_non_missing <- function(x) {
    
    idx <- which(
      !is.na(x)
    )[1]
    
    if (length(idx) == 0) {
      return(0)
    }
    
    x[idx]
  }
  
  
  a <- c(
    
    get_first_non_missing(
      y[,1]
    ),
    0,
    
    get_first_non_missing(
      y[,2]
    ),
    0,
    
    get_first_non_missing(
      y[,3]
    ),
    0
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
  # RETURN NEGATIVE LOG-LIKELIHOOD
  ##########################################################
  
  if (
    !is.finite(filt$ll)
  ) {
    
    return(1e100)
    
  }
  
  -filt$ll
}


############################################################
# FIT MODEL
############################################################

cat("\n")
cat("FITTING CTCRW MODEL\n")
cat("====================\n")

fit <- optim(
  
  par = params_start,
  
  fn = function(p) {
    neg_loglikelihood_noerror(
      p,
      aug
    )
  },
  
  method = "L-BFGS-B",
  
  control = list(
    trace = 1,
    maxit = 1000
  )
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

p_hat <- exp(
  fit$par
)


############################################################
# EXTRACT ESTIMATED PARAMETERS
############################################################

beta1_hat <- p_hat["beta1"]

beta2_hat <- p_hat["beta2"]

sigma1_hat <- p_hat["sigma1"]

sigma2_hat <- p_hat["sigma2"]


############################################################
# PRINT ESTIMATED PARAMETERS
############################################################

cat("\n")
cat("ESTIMATED PARAMETERS\n")
cat("====================\n")

cat(
  "beta1:",
  beta1_hat,
  "\n"
)

cat(
  "beta2:",
  beta2_hat,
  "\n"
)

cat(
  "sigma1:",
  sigma1_hat,
  "\n"
)

cat(
  "sigma2:",
  sigma2_hat,
  "\n"
)


############################################################
# INTERPRETABLE ESTIMATED PARAMETERS
############################################################

tau1_hat <-
  1 / beta1_hat

tau2_hat <-
  1 / beta2_hat

nu1_hat <-
  sigma1_hat /
  sqrt(
    2 * beta1_hat
  )

nu2_hat <-
  sigma2_hat /
  sqrt(
    2 * beta2_hat
  )


cat("\n")
cat("ESTIMATED INTERPRETABLE PARAMETERS\n")
cat("==================================\n")

cat(
  "Horizontal autocorrelation time tau1 (seconds): ",
  tau1_hat,
  "\n"
)

cat(
  "Horizontal autocorrelation time tau1 (hours):   ",
  tau1_hat / 3600,
  "\n"
)

cat(
  "Vertical autocorrelation time tau2 (seconds):   ",
  tau2_hat,
  "\n"
)

cat(
  "Vertical autocorrelation time tau2 (hours):     ",
  tau2_hat / 3600,
  "\n"
)

cat(
  "Horizontal RMS speed nu1 (m/s):                 ",
  nu1_hat,
  "\n"
)

cat(
  "Vertical RMS speed nu2 (m/s):                   ",
  nu2_hat,
  "\n"
)


############################################################
# TRUE VS ESTIMATED PARAMETERS
############################################################

if (
  exists("true_parameters")
) {
  
  
  ##########################################################
  # EXTRACT TRUE VALUES
  ##########################################################
  
  beta1_true <-
    true_parameters$true_value[
      true_parameters$parameter == "beta1"
    ]
  
  beta2_true <-
    true_parameters$true_value[
      true_parameters$parameter == "beta2"
    ]
  
  sigma1_true <-
    true_parameters$true_value[
      true_parameters$parameter == "sigma1"
    ]
  
  sigma2_true <-
    true_parameters$true_value[
      true_parameters$parameter == "sigma2"
    ]
  
  
  ##########################################################
  # COMPARISON TABLE
  ##########################################################
  
  comparison <- data.frame(
    
    parameter = c(
      "beta1",
      "beta2",
      "sigma1",
      "sigma2"
    ),
    
    true_value = c(
      beta1_true,
      beta2_true,
      sigma1_true,
      sigma2_true
    ),
    
    estimated_value = c(
      beta1_hat,
      beta2_hat,
      sigma1_hat,
      sigma2_hat
    )
  )
  
  
  comparison$percent_error <-
    100 *
    (
      comparison$estimated_value -
        comparison$true_value
    ) /
    comparison$true_value
  
  
  cat("\n")
  cat("TRUE VS ESTIMATED PARAMETERS\n")
  cat("============================\n")
  
  print(
    comparison,
    row.names = FALSE
  )
  
  
  ##########################################################
  # SAVE COMPARISON
  ##########################################################
  
  write.csv(
    
    comparison,
    
    "sperm_whale_CTCRW_parameter_recovery.csv",
    
    row.names = FALSE
  )
  
  
  cat("\n")
  cat(
    "Parameter recovery table saved to:\n",
    "sperm_whale_CTCRW_parameter_recovery.csv\n"
  )
  
}


############################################################
# FILTER + SMOOTHER
#
# USE ESTIMATED PARAMETERS
############################################################


############################################################
# TIME INTERVALS
############################################################

delta <- c(
  0,
  diff(aug$Time_sec)
)

delta[!is.finite(delta)] <- 0

delta <- pmax(
  delta,
  0
)


############################################################
# ZERO MEASUREMENT ERROR
############################################################

Hmat_zero <- matrix(
  0,
  N,
  3
)


############################################################
# INITIAL STATE
############################################################

get_first_non_missing <- function(x) {
  
  idx <- which(
    !is.na(x)
  )[1]
  
  if (length(idx) == 0) {
    return(0)
  }
  
  x[idx]
}


a0 <- c(
  
  get_first_non_missing(
    y[,1]
  ),
  0,
  
  get_first_non_missing(
    y[,2]
  ),
  0,
  
  get_first_non_missing(
    y[,3]
  ),
  0
)


############################################################
# INITIAL COVARIANCE
############################################################

P0 <- diag(
  c(
    100,
    10,
    100,
    10,
    100,
    10
  )
)


############################################################
# KALMAN FILTER
############################################################

cat("\n")
cat("RUNNING KALMAN FILTER...\n")

filt <- CTCRW_filter1(
  
  y = y,
  
  Hmat = Hmat_zero,
  
  beta1_vec = rep(
    beta1_hat,
    N
  ),
  
  beta2_vec = rep(
    beta2_hat,
    N
  ),
  
  s_horiz = sigma1_hat^2,
  
  s_vert = sigma2_hat^2,
  
  delta = delta,
  
  a = a0,
  
  P = P0
)


############################################################
# KALMAN SMOOTHER
############################################################

cat(
  "RUNNING KALMAN SMOOTHER...\n"
)

smooth <- CTCRW_smoother1(
  
  filter_out = filt,
  
  beta1_vec = rep(
    beta1_hat,
    N
  ),
  
  beta2_vec = rep(
    beta2_hat,
    N
  ),
  
  s_horiz = sigma1_hat^2,
  
  s_vert = sigma2_hat^2,
  
  delta = delta
)


############################################################
# CREATE SMOOTHED TRACK DATA FRAME
############################################################

smooth_track <- as.data.frame(
  smooth$a_s
)

names(smooth_track) <- c(
  "x",
  "vx",
  "y",
  "vy",
  "depth",
  "vdepth"
)

smooth_track$time <-
  aug$time


############################################################
# COMBINE OBSERVATIONS AND SMOOTHED TRACK
############################################################

plot_data <- aug %>%
  
  mutate(
    
    smooth_x =
      smooth_track$x,
    
    smooth_y =
      smooth_track$y,
    
    smooth_depth =
      smooth_track$depth
  )


############################################################
# FIRST 200 ROWS
############################################################

n_plot <- min(
  200,
  nrow(plot_data)
)

subset_plot <-
  plot_data[
    1:n_plot,
  ]


############################################################
# PLOT 1
#
# DEPTH VS TIME
############################################################

p_depth <- ggplot() +
  
  geom_point(
    
    data = subset_plot,
    
    aes(
      x = time,
      y = depth
    ),
    
    color = "blue",
    
    alpha = 0.5
  ) +
  
  geom_line(
    
    data = subset_plot,
    
    aes(
      x = time,
      y = smooth_depth
    ),
    
    color = "red",
    
    linewidth = 1
  ) +
  
  labs(
    
    title =
      paste0(
        "Depth vs Time — CTCRW Smoothed Track (First ",
        n_plot,
        " rows)"
      ),
    
    subtitle =
      "Blue = observations; red = smoothed latent trajectory",
    
    x = "Time",
    
    y = "Depth"
  ) +
  
  theme_minimal()

print(
  p_depth
)


############################################################
# PLOT 2
#
# X VS TIME
############################################################

p_x <- ggplot() +
  
  geom_point(
    
    data = subset_plot,
    
    aes(
      x = time,
      y = x
    ),
    
    color = "blue",
    
    alpha = 0.5
  ) +
  
  geom_line(
    
    data = subset_plot,
    
    aes(
      x = time,
      y = smooth_x
    ),
    
    color = "red",
    
    linewidth = 1
  ) +
  
  labs(
    
    title =
      paste0(
        "X vs Time — CTCRW Smoothed Track (First ",
        n_plot,
        " rows)"
      ),
    
    subtitle =
      "Blue = observations; red = smoothed latent trajectory",
    
    x = "Time",
    
    y = "X"
  ) +
  
  theme_minimal()

print(
  p_x
)


############################################################
# PLOT 3
#
# Y VS TIME
############################################################

p_y <- ggplot() +
  
  geom_point(
    
    data = subset_plot,
    
    aes(
      x = time,
      y = y
    ),
    
    color = "blue",
    
    alpha = 0.5
  ) +
  
  geom_line(
    
    data = subset_plot,
    
    aes(
      x = time,
      y = smooth_y
    ),
    
    color = "red",
    
    linewidth = 1
  ) +
  
  labs(
    
    title =
      paste0(
        "Y vs Time — CTCRW Smoothed Track (First ",
        n_plot,
        " rows)"
      ),
    
    subtitle =
      "Blue = observations; red = smoothed latent trajectory",
    
    x = "Time",
    
    y = "Y"
  ) +
  
  theme_minimal()

print(
  p_y
)


############################################################
# OPTIONAL:
# COMPARE SMOOTHED TRACK TO TRUE LATENT TRACK
############################################################
#
# The simulated CSV contains:
#
#   true_x
#   true_y
#   true_depth
#
# These are NOT used during fitting.
#
# They are only used here to visually assess whether
# the smoother recovered the latent trajectory.
############################################################

if (
  all(
    c(
      "true_x",
      "true_y",
      "true_depth"
    ) %in% names(aug)
  )
) {
  
  
  ##########################################################
  # TRUE VS SMOOTHED DATA
  ##########################################################
  
  comparison_track <- plot_data %>%
    
    mutate(
      
      true_x = aug$true_x,
      
      true_y = aug$true_y,
      
      true_depth = aug$true_depth
    )
  
  
  ##########################################################
  # DEPTH: TRUE VS SMOOTHED
  ##########################################################
  
  p_depth_true <- ggplot(
    comparison_track[1:n_plot, ]
  ) +
    
    geom_line(
      
      aes(
        x = time,
        y = true_depth
      ),
      
      color = "black",
      
      linewidth = 1
    ) +
    
    geom_line(
      
      aes(
        x = time,
        y = smooth_depth
      ),
      
      color = "red",
      
      linewidth = 1
    ) +
    
    labs(
      
      title =
        paste0(
          "True vs Smoothed Depth — First ",
          n_plot,
          " rows"
        ),
      
      subtitle =
        "Black = true latent trajectory; red = CTCRW smoothed trajectory",
      
      x = "Time",
      
      y = "Depth"
    ) +
    
    theme_minimal()
  
  
  print(
    p_depth_true
  )
  
  
  ##########################################################
  # X: TRUE VS SMOOTHED
  ##########################################################
  
  p_x_true <- ggplot(
    comparison_track[1:n_plot, ]
  ) +
    
    geom_line(
      
      aes(
        x = time,
        y = true_x
      ),
      
      color = "black",
      
      linewidth = 1
    ) +
    
    geom_line(
      
      aes(
        x = time,
        y = smooth_x
      ),
      
      color = "red",
      
      linewidth = 1
    ) +
    
    labs(
      
      title =
        paste0(
          "True vs Smoothed X — First ",
          n_plot,
          " rows"
        ),
      
      subtitle =
        "Black = true latent trajectory; red = CTCRW smoothed trajectory",
      
      x = "Time",
      
      y = "X"
    ) +
    
    theme_minimal()
  
  
  print(
    p_x_true
  )
  
  
  ##########################################################
  # Y: TRUE VS SMOOTHED
  ##########################################################
  
  p_y_true <- ggplot(
    comparison_track[1:n_plot, ]
  ) +
    
    geom_line(
      
      aes(
        x = time,
        y = true_y
      ),
      
      color = "black",
      
      linewidth = 1
    ) +
    
    geom_line(
      
      aes(
        x = time,
        y = smooth_y
      ),
      
      color = "red",
      
      linewidth = 1
    ) +
    
    labs(
      
      title =
        paste0(
          "True vs Smoothed Y — First ",
          n_plot,
          " rows"
        ),
      
      subtitle =
        "Black = true latent trajectory; red = CTCRW smoothed trajectory",
      
      x = "Time",
      
      y = "Y"
    ) +
    
    theme_minimal()
  
  
  print(
    p_y_true
  )
  
  
  ##########################################################
  # SAVE TRUE VS SMOOTHED TRACK
  ##########################################################
  
  write.csv(
    
    comparison_track,
    
    "sperm_whale_CTCRW_true_vs_smoothed.csv",
    
    row.names = FALSE
  )
  
}


############################################################
# SAVE SMOOTHED TRACK
############################################################

write.csv(
  
  smooth_track,
  
  "sperm_whale_CTCRW_smoothed_track.csv",
  
  row.names = FALSE
)


############################################################
# FINAL SUMMARY
############################################################

cat("\n")
cat("====================================================\n")
cat("SIMULATION-RECOVERY TEST COMPLETE\n")
cat("====================================================\n")

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

cat("\n")

cat(
  "Estimated beta1: ",
  beta1_hat,
  "\n"
)

cat(
  "Estimated beta2: ",
  beta2_hat,
  "\n"
)

cat(
  "Estimated sigma1: ",
  sigma1_hat,
  "\n"
)

cat(
  "Estimated sigma2: ",
  sigma2_hat,
  "\n"
)

cat("\n")

cat(
  "Smoothed track saved to:\n",
  "sperm_whale_CTCRW_smoothed_track.csv\n"
)

cat("\n")

cat("DONE.\n")










