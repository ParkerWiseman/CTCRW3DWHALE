


############################################################
# CTCRW 3D SPERM WHALE MODEL — ABSOLUTELY NO MEASUREMENT ERROR
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

tau_horiz_hours <- 8
tau_vert_hours  <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours  * 3600

beta1_start <- 1 / tau_horiz_sec
beta2_start <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

sigma1_start <- v_rms_horiz * sqrt(2 * beta1_start)
sigma2_start <- v_rms_vert  * sqrt(2 * beta2_start)

params_start <- c(
  beta1  = log(beta1_start),
  beta2  = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

print(params_start)

############################################################
# MODIFY NEGATIVE LOG-LIKELIHOOD TO USE ZERO ERROR
############################################################

neg_loglikelihood_noerror <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  
  delta_raw <- c(NA, diff(time_sec))
  delta_raw[1] <- delta_raw[2]
  
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### ZERO MEASUREMENT ERROR
  Hmat <- matrix(0, nrow(data_aug), 3)
  
  get_first_non_missing <- function(col) {
    idx <- which(!is.na(col))[1]
    if (is.na(idx)) 0 else col[idx]
  }
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  P <- diag(6) * 1e6
  
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
# FIT MODEL
############################################################

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood_noerror(p, aug),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta1_hat
tau2_hat <- 1 / beta2_hat

nu1_hat <- sigma1_hat / sqrt(2 * beta1_hat)
nu2_hat <- sigma2_hat / sqrt(2 * beta2_hat)

cat("Horizontal autocorrelation time tau1 (s):", tau1_hat, "\n")
cat("Vertical autocorrelation time tau2 (s):", tau2_hat, "\n")
cat("Horizontal RMS speed nu1 (m/s):", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):", nu2_hat, "\n")

############################################################
# FILTER AND SMOOTHER
############################################################

Hmat_zero <- matrix(0, N, 3)

delta_raw <- c(NA, diff(aug$Time_sec))
delta_raw[1] <- delta_raw[2]
delta <- pmax(delta_raw, 1)
delta[!is.finite(delta)] <- 1

a0 <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
P0 <- diag(c(100, 10, 100, 10, 100, 10))


filt <- CTCRW_filter1(
  y         = y,
  Hmat      = Hmat_zero,
  beta1_vec = rep(beta1_hat, N),
  beta2_vec = rep(beta2_hat, N),
  s_horiz   = sigma1_hat^2,
  s_vert    = sigma2_hat^2,
  delta     = delta,
  a         = a0,
  P         = P0
)

smooth <- CTCRW_smoother1(
  filter_out = filt,
  beta1_vec  = rep(beta1_hat, N),
  beta2_vec  = rep(beta2_hat, N),
  s_horiz    = sigma1_hat^2,
  s_vert     = sigma2_hat^2,
  delta      = delta
)

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time



############################################################
# PLOTS FOR FIRST 200 ROWS — WITH SMOOTHED VALUES
############################################################

subset_obs    <- aug[1:200, ]
subset_smooth <- smooth_track[1:200, ]

## 1. Depth vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = depth),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = depth),
            color = "red", linewidth = 1) +
  labs(title = "Depth vs Time (First 200 rows)",
       x = "Time", y = "Depth") +
  theme_minimal()

## 2. X vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = x),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = x),
            color = "red", linewidth = 1) +
  labs(title = "X vs Time (First 200 rows)",
       x = "Time", y = "X") +
  theme_minimal()

## 3. Y vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = y),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = y),
            color = "red", linewidth = 1) +
  labs(title = "Y vs Time (First 200 rows)",
       x = "Time", y = "Y") +
  theme_minimal()

































############################################################
# CTCRW 3D SPERM WHALE MODEL — BASICALLY NO MEASUREMENT ERROR & estimating all parameters
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

tau_horiz_hours <- 8
tau_vert_hours  <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours  * 3600

beta1_start <- 1 / tau_horiz_sec
beta2_start <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

sigma1_start <- v_rms_horiz * sqrt(2 * beta1_start)
sigma2_start <- v_rms_vert  * sqrt(2 * beta2_start)

params_start <- c(
  beta1  = log(beta1_start),
  beta2  = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

print(params_start)

############################################################
# MODIFY NEGATIVE LOG-LIKELIHOOD TO USE ZERO ERROR
############################################################

neg_loglikelihood_noerror <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  
  delta_raw <- c(NA, diff(time_sec))
  delta_raw[1] <- delta_raw[2]
  
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### ZERO MEASUREMENT ERROR
  Hmat <- matrix(0, nrow(data_aug), 3)
  
  get_first_non_missing <- function(col) {
    idx <- which(!is.na(col))[1]
    if (is.na(idx)) 0 else col[idx]
  }
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  P <- diag(6) * 1e6
  
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
# FIT MODEL
############################################################

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood_noerror(p, aug),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta1_hat
tau2_hat <- 1 / beta2_hat

nu1_hat <- sigma1_hat / sqrt(2 * beta1_hat)
nu2_hat <- sigma2_hat / sqrt(2 * beta2_hat)

cat("Horizontal autocorrelation time tau1 (s):", tau1_hat, "\n")
cat("Vertical autocorrelation time tau2 (s):", tau2_hat, "\n")
cat("Horizontal RMS speed nu1 (m/s):", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):", nu2_hat, "\n")

############################################################
# FILTER AND SMOOTHER
############################################################

Hmat_zero <- matrix(1e-6, N, 3)

delta_raw <- c(NA, diff(aug$Time_sec))
delta_raw[1] <- delta_raw[2]
delta <- pmax(delta_raw, 1)
delta[!is.finite(delta)] <- 1

a0 <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
P0 <- diag(c(100, 10, 100, 10, 100, 10))


filt <- CTCRW_filter1(
  y         = y,
  Hmat      = Hmat_zero,
  beta1_vec = rep(beta1_hat, N),
  beta2_vec = rep(beta2_hat, N),
  s_horiz   = sigma1_hat^2,
  s_vert    = sigma2_hat^2,
  delta     = delta,
  a         = a0,
  P         = P0
)

smooth <- CTCRW_smoother1(
  filter_out = filt,
  beta1_vec  = rep(beta1_hat, N),
  beta2_vec  = rep(beta2_hat, N),
  s_horiz    = sigma1_hat^2,
  s_vert     = sigma2_hat^2,
  delta      = delta
)

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time




############################################################
# PLOTS FOR FIRST 1000 ROWS — WITH SMOOTHED VALUES
############################################################

subset_obs    <- aug[1:1000, ]
subset_smooth <- smooth_track[1:1000, ]

# Remove rows where smoother produced NA values
subset_smooth <- subset_smooth %>% 
  filter(!is.na(x), !is.na(y), !is.na(depth))

## 1. Depth vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = depth),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = depth),
            color = "red", linewidth = 1) +
  labs(title = "Depth vs Time (First 1000 rows)",
       x = "Time", y = "Depth") +
  theme_minimal()

## 2. X vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = x),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = x),
            color = "red", linewidth = 1) +
  labs(title = "X vs Time (First 1000 rows)",
       x = "Time", y = "X") +
  theme_minimal()

## 3. Y vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = y),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = y),
            color = "red", linewidth = 1) +
  labs(title = "Y vs Time (First 1000 rows)",
       x = "Time", y = "Y") +
  theme_minimal()















############################################################
# CTCRW 3D SPERM WHALE MODEL — BASICALLY NO MEASUREMENT ERROR & FIXING BETAS
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
# MICHELOT & BLACKWELL (2020) VALUES
############################################################

tau_horiz_hours <- 8
tau_vert_hours  <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours  * 3600

beta1_fixed <- 1 / tau_horiz_sec
beta2_fixed <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

sigma1_start <- v_rms_horiz * sqrt(2 * beta1_fixed)
sigma2_start <- v_rms_vert  * sqrt(2 * beta2_fixed)

params_start <- c(
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

print(params_start)

############################################################
# NEGATIVE LOG-LIKELIHOOD — β FIXED, σ ESTIMATED
############################################################

neg_loglikelihood_fixbeta <- function(params, data_aug) {
  
  # FIXED BETAS
  beta1 <- beta1_fixed
  beta2 <- beta2_fixed
  
  # ESTIMATED SIGMAS
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  
  delta_raw <- c(NA, diff(time_sec))
  delta_raw[1] <- delta_raw[2]
  
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  # ZERO MEASUREMENT ERROR
  Hmat <- matrix(1e-6, nrow(data_aug), 3)
  
  get_first_non_missing <- function(col) {
    idx <- which(!is.na(col))[1]
    if (is.na(idx)) 0 else col[idx]
  }
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  P <- diag(6) * 1e6
  
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
# FIT MODEL — ONLY σ₁ AND σ₂ ARE ESTIMATED
############################################################

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood_fixbeta(p, aug),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)

sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

beta1_hat <- beta1_fixed
beta2_hat <- beta2_fixed

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta1_hat
tau2_hat <- 1 / beta2_hat

nu1_hat <- sigma1_hat / sqrt(2 * beta1_hat)
nu2_hat <- sigma2_hat / sqrt(2 * beta2_hat)

cat("Horizontal autocorrelation time tau1 (s):", tau1_hat, "\n")
cat("Vertical autocorrelation time tau2 (s):", tau2_hat, "\n")
cat("Horizontal RMS speed nu1 (m/s):", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):", nu2_hat, "\n")

############################################################
# FILTER AND SMOOTHER
############################################################

Hmat_zero <- matrix(1e-6, N, 3)

delta_raw <- c(NA, diff(aug$Time_sec))
delta_raw[1] <- delta_raw[2]
delta <- pmax(delta_raw, 1)
delta[!is.finite(delta)] <- 1

a0 <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
P0 <- diag(c(100, 10, 100, 10, 100, 10))

filt <- CTCRW_filter1(
  y         = y,
  Hmat      = Hmat_zero,
  beta1_vec = rep(beta1_hat, N),
  beta2_vec = rep(beta2_hat, N),
  s_horiz   = sigma1_hat^2,
  s_vert    = sigma2_hat^2,
  delta     = delta,
  a         = a0,
  P         = P0
)

smooth <- CTCRW_smoother1(
  filter_out = filt,
  beta1_vec  = rep(beta1_hat, N),
  beta2_vec  = rep(beta2_hat, N),
  s_horiz    = sigma1_hat^2,
  s_vert     = sigma2_hat^2,
  delta      = delta
)

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time

############################################################
# PLOTS FOR FIRST 1000 ROWS
############################################################

subset_obs    <- aug[1:1000, ]
subset_smooth <- smooth_track[1:1000, ]

subset_smooth <- subset_smooth %>% 
  filter(!is.na(x), !is.na(y), !is.na(depth))

ggplot() +
  geom_point(data = subset_obs, aes(time, depth), color="blue", alpha=0.5) +
  geom_line(data = subset_smooth, aes(time, depth), color="red", linewidth=1) +
  labs(title="Depth vs Time (First 1000 rows)", x="Time", y="Depth") +
  theme_minimal()

ggplot() +
  geom_point(data = subset_obs, aes(time, x), color="blue", alpha=0.5) +
  geom_line(data = subset_smooth, aes(time, x), color="red", linewidth=1) +
  labs(title="X vs Time (First 1000 rows)", x="Time", y="X") +
  theme_minimal()

ggplot() +
  geom_point(data = subset_obs, aes(time, y), color="blue", alpha=0.5) +
  geom_line(data = subset_smooth, aes(time, y), color="red", linewidth=1) +
  labs(title="Y vs Time (First 1000 rows)", x="Time", y="Y") +
  theme_minimal()















# good one

############################################################
# CTCRW 3D SPERM WHALE MODEL — tiny measurement error & estimating all parameters
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
    Time_sec   = as.numeric(difftime(time, min(time), units = "secs")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

############################################################
# INITIAL VALUES FROM MICHELOT & BLACKWELL (2020)
############################################################

tau_horiz_hours <- 8
tau_vert_hours  <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours  * 3600

beta1_start <- 1 / tau_horiz_sec
beta2_start <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

sigma1_start <- v_rms_horiz * sqrt(2 * beta1_start)
sigma2_start <- v_rms_vert  * sqrt(2 * beta2_start)

params_start <- c(
  beta1  = log(beta1_start),
  beta2  = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

print(params_start)

############################################################
# NEGATIVE LOG-LIKELIHOOD WITH TINY MEASUREMENT ERROR
############################################################

neg_loglikelihood_noerror <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  
  delta_raw <- c(NA, diff(time_sec))
  delta_raw[1] <- delta_raw[2]
  
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### TINY MEASUREMENT ERROR (stabilizes filter)
  Hmat <- matrix(0, nrow(data_aug), 3)
  Hmat[,1] <- 5^2
  Hmat[,2] <- 5^2
  Hmat[,3] <- 10^2
  
  get_first_non_missing <- function(col) {
    idx <- which(!is.na(col))[1]
    if (is.na(idx)) 0 else col[idx]
  }
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  P <- diag(6) * 1e6
  
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
# FIT MODEL (estimate beta1, beta2, sigma1, sigma2)
############################################################

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood_noerror(p, aug),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta1_hat
tau2_hat <- 1 / beta2_hat

nu1_hat <- sigma1_hat / sqrt(2 * beta1_hat)
nu2_hat <- sigma2_hat / sqrt(2 * beta2_hat)

cat("Horizontal autocorrelation time tau1 (s):", tau1_hat, "\n")
cat("Vertical autocorrelation time tau2 (s):", tau2_hat, "\n")
cat("Horizontal RMS speed nu1 (m/s):", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):", nu2_hat, "\n")

############################################################
# FILTER AND SMOOTHER WITH SAME TINY MEASUREMENT ERROR
############################################################

Hmat_zero <- matrix(0, N, 3)
Hmat_zero[,1] <- 5^2
Hmat_zero[,2] <- 5^2
Hmat_zero[,3] <- 10^2


delta_raw <- c(NA, diff(aug$Time_sec))
delta_raw[1] <- delta_raw[2]
delta <- pmax(delta_raw, 1)
delta[!is.finite(delta)] <- 1

a0 <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
P0 <- diag(c(100, 10, 100, 10, 100, 10))

filt <- CTCRW_filter1(
  y         = y,
  Hmat      = Hmat_zero,
  beta1_vec = rep(beta1_hat, N),
  beta2_vec = rep(beta2_hat, N),
  s_horiz   = sigma1_hat^2,
  s_vert    = sigma2_hat^2,
  delta     = delta,
  a         = a0,
  P         = P0
)

smooth <- CTCRW_smoother1(
  filter_out = filt,
  beta1_vec  = rep(beta1_hat, N),
  beta2_vec  = rep(beta2_hat, N),
  s_horiz    = sigma1_hat^2,
  s_vert     = sigma2_hat^2,
  delta      = delta
)

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time

############################################################
# PLOTS FOR FIRST 1000 ROWS — WITH SMOOTHED VALUES
############################################################

subset_obs    <- aug[1:1000, ]
subset_smooth <- smooth_track[1:1000, ]

subset_smooth <- subset_smooth %>% 
  filter(!is.na(x), !is.na(y), !is.na(depth))

## 1. Depth vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = depth),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = depth),
            color = "red", linewidth = 1) +
  labs(title = "Depth vs Time (First 1000 rows)",
       x = "Time", y = "Depth") +
  theme_minimal()

## 2. X vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = x),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = x),
            color = "red", linewidth = 1) +
  labs(title = "X vs Time (First 1000 rows)",
       x = "Time", y = "X") +
  theme_minimal()

## 3. Y vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = y),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = y),
            color = "red", linewidth = 1) +
  labs(title = "Y vs Time (First 1000 rows)",
       x = "Time", y = "Y") +
  theme_minimal()

















############################################################
# CTCRW 3D SPERM WHALE MODEL — tiny measurement error & estimating all parameters
# & changing the way time intervals are computed
# (using delta_raw <- c(diff(time_sec), 1) instead of delta_raw <- c(NA, diff(time_sec)))
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
    Time_sec   = as.numeric(difftime(time, min(time), units = "secs")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

############################################################
# INITIAL VALUES FROM MICHELOT & BLACKWELL (2020)
############################################################

tau_horiz_hours <- 8
tau_vert_hours  <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours  * 3600

beta1_start <- 1 / tau_horiz_sec
beta2_start <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

sigma1_start <- v_rms_horiz * sqrt(2 * beta1_start)
sigma2_start <- v_rms_vert  * sqrt(2 * beta2_start)

params_start <- c(
  beta1  = log(beta1_start),
  beta2  = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

print(params_start)

############################################################
# NEGATIVE LOG-LIKELIHOOD WITH TINY MEASUREMENT ERROR
############################################################

neg_loglikelihood_noerror <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  
  delta_raw <- c(diff(time_sec), 1)
  delta_raw[1] <- delta_raw[2]
  
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### TINY MEASUREMENT ERROR (stabilizes filter)
  Hmat <- matrix(0, nrow(data_aug), 3)
  Hmat[,1] <- 5^2
  Hmat[,2] <- 5^2
  Hmat[,3] <- 10^2
  
  get_first_non_missing <- function(col) {
    idx <- which(!is.na(col))[1]
    if (is.na(idx)) 0 else col[idx]
  }
  
  a <- c(
    get_first_non_missing(y[,1]), 0,
    get_first_non_missing(y[,2]), 0,
    get_first_non_missing(y[,3]), 0
  )
  
  P <- diag(6) * 1e6
  
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
# FIT MODEL (estimate beta1, beta2, sigma1, sigma2)
############################################################

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood_noerror(p, aug),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

############################################################
# INTERPRETABLE PARAMETERS
############################################################

tau1_hat <- 1 / beta1_hat
tau2_hat <- 1 / beta2_hat

nu1_hat <- sigma1_hat / sqrt(2 * beta1_hat)
nu2_hat <- sigma2_hat / sqrt(2 * beta2_hat)

cat("Horizontal autocorrelation time tau1 (s):", tau1_hat, "\n")
cat("Vertical autocorrelation time tau2 (s):", tau2_hat, "\n")
cat("Horizontal RMS speed nu1 (m/s):", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):", nu2_hat, "\n")

############################################################
# FILTER AND SMOOTHER WITH SAME TINY MEASUREMENT ERROR
############################################################

Hmat_zero <- matrix(0, N, 3)
Hmat_zero[,1] <- 5^2
Hmat_zero[,2] <- 5^2
Hmat_zero[,3] <- 10^2


delta_raw <- c(NA, diff(aug$Time_sec))
delta_raw[1] <- delta_raw[2]
delta <- pmax(delta_raw, 1)
delta[!is.finite(delta)] <- 1

a0 <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
P0 <- diag(c(100, 10, 100, 10, 100, 10))

filt <- CTCRW_filter1(
  y         = y,
  Hmat      = Hmat_zero,
  beta1_vec = rep(beta1_hat, N),
  beta2_vec = rep(beta2_hat, N),
  s_horiz   = sigma1_hat^2,
  s_vert    = sigma2_hat^2,
  delta     = delta,
  a         = a0,
  P         = P0
)

smooth <- CTCRW_smoother1(
  filter_out = filt,
  beta1_vec  = rep(beta1_hat, N),
  beta2_vec  = rep(beta2_hat, N),
  s_horiz    = sigma1_hat^2,
  s_vert     = sigma2_hat^2,
  delta      = delta
)

smooth_track <- as.data.frame(smooth$a_s)
names(smooth_track) <- c("x","vx","y","vy","depth","vdepth")
smooth_track$time <- aug$time

############################################################
# PLOTS FOR FIRST 1000 ROWS — WITH SMOOTHED VALUES
############################################################

subset_obs    <- aug[1:1000, ]
subset_smooth <- smooth_track[1:1000, ]

subset_smooth <- subset_smooth %>% 
  filter(!is.na(x), !is.na(y), !is.na(depth))

## 1. Depth vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = depth),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = depth),
            color = "red", linewidth = 1) +
  labs(title = "Depth vs Time (First 1000 rows)",
       x = "Time", y = "Depth") +
  theme_minimal()

## 2. X vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = x),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = x),
            color = "red", linewidth = 1) +
  labs(title = "X vs Time (First 1000 rows)",
       x = "Time", y = "X") +
  theme_minimal()

## 3. Y vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = y),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = y),
            color = "red", linewidth = 1) +
  labs(title = "Y vs Time (First 1000 rows)",
       x = "Time", y = "Y") +
  theme_minimal()






##########################################

# September 16, 2026



############################################################
# CTCRW 3D SPERM WHALE MODEL — tiny measurement error & estimating all parameters
# & changing the way time intervals are computed
# (using delta_raw <- c(diff(time_sec), 1) instead of delta_raw <- c(NA, diff(time_sec)))
############################################################



############################################################
# CTCRW 3D SPERM WHALE MODEL — linear HDOP error + estimating all parameters
# using delta_raw <- c(diff(time_sec), 1)
# Instead of using neg_loglikelihood_noerror(), we use
# neg_loglikelihood(..., error_model = "linearerror"), which internally calls build_Hmat_LinearError().
# Also, in the filter/smooother H matrix portion of this code we also instead
# use Hmat <- build_Hmat_LinearError(...)
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
    Time_sec   = as.numeric(difftime(time, min(time), units = "secs")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

############################################################
# INITIAL VALUES FROM MICHELOT & BLACKWELL (2020)
############################################################

tau_horiz_hours <- 8
tau_vert_hours  <- 3

tau_horiz_sec <- tau_horiz_hours * 3600
tau_vert_sec  <- tau_vert_hours  * 3600

beta1_start <- 1 / tau_horiz_sec
beta2_start <- 1 / tau_vert_sec

v_rms_horiz <- 1.0
v_rms_vert  <- 0.5

sigma1_start <- v_rms_horiz * sqrt(2 * beta1_start)
sigma2_start <- v_rms_vert  * sqrt(2 * beta2_start)

params_start <- c(
  beta1  = log(beta1_start),
  beta2  = log(beta2_start),
  sigma1 = log(sigma1_start),
  sigma2 = log(sigma2_start)
)

print(params_start)

############################################################
# FIT MODEL — using linear HDOP error model
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

############################################################
# INTERPRETABLE PARAMETERS — tau and nu
############################################################

tau1_hat <- 1 / beta1_hat
tau2_hat <- 1 / beta2_hat

nu1_hat <- sigma1_hat / sqrt(2 * beta1_hat)
nu2_hat <- sigma2_hat / sqrt(2 * beta2_hat)

cat("Horizontal autocorrelation time tau1 (s):", tau1_hat, "\n")
cat("Vertical autocorrelation time tau2 (s):", tau2_hat, "\n")
cat("Horizontal RMS speed nu1 (m/s):", nu1_hat, "\n")
cat("Vertical RMS speed nu2 (m/s):", nu2_hat, "\n")

s_horiz_hat <- sigma1_hat^2
s_vert_hat  <- sigma2_hat^2

############################################################
# BUILD H MATRIX — linear HDOP model
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

delta_raw <- c(diff(aug$Time_sec), 1)
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

P0 <- diag(c(100, 10, 100, 10, 100, 10))

############################################################
# FILTER — using linear HDOP Hmat
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
# SMOOTHER — using linear HDOP Hmat
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
# PLOTS FOR FIRST 1000 ROWS — WITH SMOOTHED VALUES
############################################################

subset_obs    <- aug[1:1000, ]
subset_smooth <- smooth_track[1:1000, ]

subset_smooth <- subset_smooth %>% 
  filter(!is.na(x), !is.na(y), !is.na(depth))

## 1. Depth vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = depth),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = depth),
            color = "red", linewidth = 1) +
  labs(title = "Depth vs Time (First 1000 rows)",
       x = "Time", y = "Depth") +
  theme_minimal()

## 2. X vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = x),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = x),
            color = "red", linewidth = 1) +
  labs(title = "X vs Time (First 1000 rows)",
       x = "Time", y = "X") +
  theme_minimal()

## 3. Y vs Time
ggplot() +
  geom_point(data = subset_obs,
             aes(x = time, y = y),
             color = "blue", alpha = 0.5) +
  geom_line(data = subset_smooth,
            aes(x = time, y = y),
            color = "red", linewidth = 1) +
  labs(title = "Y vs Time (First 1000 rows)",
       x = "Time", y = "Y") +
  theme_minimal()















