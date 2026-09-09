


# Using LinearHDOP code to analyze the sperm whale datatset, obtain parameter estimates,
# and obtain a smoothed movement track.



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
source("simulate_CTCRW_3D.R") # still useful for testing

############################################################
# LOAD REAL 3D SPERM WHALE DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")   # <-- FIXED

whale$time <- ymd_hms(whale$time)

aug <- whale %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)
dt <- median(diff(aug$Time))   # time step in days

############################################################
# PARAMETER ESTIMATION USING LINEAR HDOP MODEL
############################################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
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

p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

s_horiz_hat <- sigma1_hat^2
s_vert_hat  <- sigma2_hat^2

############################################################
# BUILD H MATRIX FOR WHALE DATA (LINEAR HDOP)
############################################################

var0_xy  <- 0
var1_xy  <- 0.1
sd_depth <- 10

Hmat <- build_Hmat_LinearError2(aug, var0_xy, var1_xy, sd_depth)

############################################################
# INITIAL STATE
############################################################

a0 <- c(
  ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
  ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
  ifelse(is.na(y[1,3]), 0, y[1,3]), 0
)
P0 <- diag(6) * 1e2

delta <- rep(dt, N)

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
# PLOT 2D (XY)
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = y), color = "red", size = 1) +
  labs(title = "Observed vs Smoothed Whale Track (XY)",
       x = "X", y = "Y") +
  theme_minimal()

############################################################
# PLOT DEPTH VS X
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = depth), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = depth), color = "red", size = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Depth vs X)",
       x = "X", y = "Depth") +
  theme_minimal()

############################################################
# 3D INTERACTIVE PLOT
############################################################

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




#################################################

# August 12, 2026




library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

source("matrices.R")          # makeT(), makeQ()
source("CTCRW_filter.R")      # CTCRW_filter1()
source("CTCRW_smoother.R")    # CTCRW_smoother1()
source("neg_loglikelihood.R") # unified neg_loglikelihood()
source("simulate_CTCRW_3D.R") # still useful for testing

############################################################
# LOAD REAL 3D SPERM WHALE DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")

whale$time <- ymd_hms(whale$time)

aug <- whale %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

############################################################
# PARAMETER ESTIMATION USING UNIFIED LINEAR HDOP MODEL
############################################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

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
# BUILD H MATRIX FOR WHALE DATA (LINEAR HDOP)
############################################################

var0_xy  <- 0
var1_xy  <- 0.1
sd_depth <- 10

Hmat <- build_Hmat_LinearError2(aug, var0_xy, var1_xy, sd_depth)

############################################################
# DELTA (MUST MATCH LIKELIHOOD)
############################################################

delta <- c(NA, diff(aug$Time))
delta[1] <- delta[2]

############################################################
# INITIAL STATE
############################################################

a0 <- c(
  ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
  ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
  ifelse(is.na(y[1,3]), 0, y[1,3]), 0
)
P0 <- diag(6) * 1e2

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
# PLOT 2D (XY)
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = y), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (XY)",
       x = "X", y = "Y") +
  theme_minimal()

############################################################
# PLOT DEPTH VS X
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = depth), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = depth), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Depth vs X)",
       x = "X", y = "Depth") +
  theme_minimal()

############################################################
# 3D INTERACTIVE PLOT
############################################################

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











###############################

# August 23, 2026



############################################################
# CTCRW 3D SPERM WHALE MODEL — CLEAN FULL SCRIPT
############################################################

library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)

set.seed(123)

############################################################
# LOAD CLEAN FUNCTION FILES
############################################################

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")
source("simulate_CTCRW_3D.R")

############################################################
# LOAD REAL SPERM WHALE DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")
whale$time <- ymd_hms(whale$time)

aug <- whale %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

############################################################
# PARAMETER ESTIMATION USING ROBUST LINEAR HDOP MODEL
############################################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood(p, aug, error_model = "linearerror"),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

#fit <- optim(
#  par      = params_start,
#  fn       = function(p) neg_loglikelihood(p, aug, error_model = "linearerror"),
#  method   = "L-BFGS-B",
#  lower    = c(log(0.01), log(0.01), log(0.1), log(0.1)),
#  upper    = c(log(100),  log(100),  log(5000), log(5000)),
#  control  = list(trace = 1, maxit = 1000)
#)


p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

s_horiz_hat <- sigma1_hat^2
s_vert_hat  <- sigma2_hat^2

############################################################
# BUILD H MATRIX (ROBUST LINEAR HDOP)
############################################################

var0_xy  <- 0
var1_xy  <- 0.1
sd_depth <- 10

Hmat <- build_Hmat_LinearError(
  aug,
  var0_xy = 25,     # baseline variance (5 m SD)
  var1_xy = 5,      # HDOP scaling
  sd_depth = 10
)
### Hmat <- build_Hmat_LinearError(aug, var0_xy, var1_xy, sd_depth)

############################################################
# DELTA (MATCHES LIKELIHOOD)
############################################################

delta <- c(NA, diff(aug$Time))
delta[1] <- delta[2]
delta <- pmax(delta, 1e-6)
delta[!is.finite(delta)] <- 1e-6

############################################################
# INITIAL STATE
############################################################

a0 <- c(
  ifelse(is.na(y[1,1]), 0, y[1,1]), 0,
  ifelse(is.na(y[1,2]), 0, y[1,2]), 0,
  ifelse(is.na(y[1,3]), 0, y[1,3]), 0
)

P0 <- diag(6) * 1e2

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
# PLOT 2D XY
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = y), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (XY)",
       x = "X", y = "Y") +
  theme_minimal()

############################################################
# PLOT DEPTH VS X
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = depth), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = depth), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Depth vs X)",
       x = "X", y = "Depth") +
  theme_minimal()

############################################################
# PLOT X VS TIME
############################################################

ggplot() +
  geom_point(data = aug, aes(x = time, y = x), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = time, y = x), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (X vs Time)",
       x = "Time", y = "X") +
  theme_minimal()

############################################################
# PLOT Y VS TIME
############################################################

ggplot() +
  geom_point(data = aug, aes(x = time, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = time, y = y), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Y vs Time)",
       x = "Time", y = "Y") +
  theme_minimal()

############################################################
# 3D PLOT
############################################################

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






##################################





# September 4, 2026

############################################################
# CTCRW 3D SPERM WHALE MODEL — CLEAN FULL SCRIPT
############################################################

library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)

set.seed(123)

############################################################
# LOAD CLEAN FUNCTION FILES
############################################################

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")
source("simulate_CTCRW_3D.R")

############################################################
# LOAD REAL SPERM WHALE DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")
whale$time <- ymd_hms(whale$time)

aug <- whale %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

############################################################
# PARAMETER ESTIMATION USING ROBUST LINEAR HDOP MODEL
############################################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood(p, aug, error_model = "linearerror"),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

#fit <- optim(
#  par      = params_start,
#  fn       = function(p) neg_loglikelihood(p, aug, error_model = "linearerror"),
#  method   = "L-BFGS-B",
#  lower    = c(log(0.01), log(0.01), log(0.1), log(0.1)),
#  upper    = c(log(100),  log(100),  log(5000), log(5000)),
#  control  = list(trace = 1, maxit = 1000)
#)


p_hat <- exp(fit$par)
print(p_hat)

beta1_hat  <- p_hat["beta1"]
beta2_hat  <- p_hat["beta2"]
sigma1_hat <- p_hat["sigma1"]
sigma2_hat <- p_hat["sigma2"]

s_horiz_hat <- sigma1_hat^2
s_vert_hat  <- sigma2_hat^2

############################################################
# BUILD H MATRIX (ROBUST LINEAR HDOP)
############################################################

var0_xy  <- 0
var1_xy  <- 0.1
sd_depth <- 10

Hmat <- build_Hmat_LinearError(
  aug,
  var0_xy = 25,     # baseline variance (5 m SD)
  var1_xy = 5,      # HDOP scaling
  sd_depth = 10
)
### Hmat <- build_Hmat_LinearError(aug, var0_xy, var1_xy, sd_depth)

############################################################
# DELTA (MATCHES LIKELIHOOD)
############################################################

delta <- c(NA, diff(aug$Time))
delta[1] <- delta[2]
delta <- pmax(delta, 1e-6)
delta[!is.finite(delta)] <- 1e-6

############################################################
# INITIAL STATE
############################################################


# Find first non‑missing x, y, depth
x0 <- y[which(!is.na(y[,1]))[1], 1]
y0 <- y[which(!is.na(y[,2]))[1], 2]
d0 <- y[which(!is.na(y[,3]))[1], 3]

# Initial state vector
a0 <- c(
  x0, 0,
  y0, 0,
  d0, 0
)

P0 <- diag(6) * 1e2

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
# PLOT 2D XY
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = y), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (XY)",
       x = "X", y = "Y") +
  theme_minimal()

############################################################
# PLOT DEPTH VS X
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = depth), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = depth), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Depth vs X)",
       x = "X", y = "Depth") +
  theme_minimal()

############################################################
# PLOT X VS TIME
############################################################

ggplot() +
  geom_point(data = aug, aes(x = time, y = x), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = time, y = x), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (X vs Time)",
       x = "Time", y = "X") +
  theme_minimal()

############################################################
# PLOT Y VS TIME
############################################################

ggplot() +
  geom_point(data = aug, aes(x = time, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = time, y = y), color = "red", linewidth = 1) +
  labs(title = "Observed vs Smoothed Whale Track (Y vs Time)",
       x = "Time", y = "Y") +
  theme_minimal()

############################################################
# 3D PLOT
############################################################

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







###################
###################

# include neg_loglikelihood code and sperm whale code in each entry.

# september 5, 2026


############################################################
# neg_loglikelihood.R — Stable version with delta in seconds
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
  #### H MATRIX
  #### ----------------------------------------------------
  
  var0_xy  <- 100
  var1_xy  <- 2
  sd_depth <- 25
  
  Hmat <- build_Hmat_LinearError(data_aug, var0_xy, var1_xy, sd_depth)
  
  #### ----------------------------------------------------
  #### INITIAL STATE
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
# CTCRW 3D SPERM WHALE MODEL — delta in seconds
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

whale <- read.csv("sperm_whale_processed.csv")
whale$time <- ymd_hms(whale$time)

aug <- whale %>%
  mutate(
    Time_sec = as.numeric(difftime(time, min(time), units = "secs")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])
N <- nrow(aug)

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

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

var0_xy  <- 100
var1_xy  <- 2
sd_depth <- 25

Hmat <- build_Hmat_LinearError(
  aug,
  var0_xy = var0_xy,
  var1_xy = var1_xy,
  sd_depth = sd_depth
)

delta_raw <- c(NA, diff(aug$Time_sec))
delta_raw[1] <- delta_raw[2]

delta <- pmax(delta_raw, 1)
delta[!is.finite(delta)] <- 1

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

# plots omitted for brevity














##########################

# September 6, 2026

# using Michelot & Blackwell formulas

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









###########################

# September 8, 2026


############################################################
# neg_loglikelihood.R — CTCRW likelihood with delta in seconds
# Supports: "noerror", "constanterror", "linearerror"
############################################################

#### ----------------------------------------------------
#### H-MATRIX BUILDERS
#### ----------------------------------------------------

build_Hmat_NoError <- function(data_aug) {
  N <- nrow(data_aug)
  eps <- 1e-6
  Hmat <- matrix(eps, N, 3)
  Hmat
}

build_Hmat_ConstantError <- function(data_aug, sd_xy, sd_depth) {
  N <- nrow(data_aug)
  Hmat <- matrix(sd_xy^2, N, 3)
  Hmat[,3] <- sd_depth^2
  Hmat
}

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

#### ----------------------------------------------------
#### NEGATIVE LOG-LIKELIHOOD
#### ----------------------------------------------------

neg_loglikelihood <- function(params, data_aug,
                              error_model = c("noerror","constanterror","linearerror")) {
  
  error_model <- match.arg(error_model)
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  #### DELTA IN SECONDS
  time_sec <- as.numeric(difftime(data_aug$time,
                                  min(data_aug$time),
                                  units = "secs"))
  delta_raw <- c(NA, diff(time_sec))
  delta_raw[1] <- delta_raw[2]
  delta <- pmax(delta_raw, 1)
  delta[!is.finite(delta)] <- 1
  
  #### PROCESS NOISE
  s_horiz <- sigma1^2
  s_vert  <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  #### H MATRIX SWITCH
  if (error_model == "noerror") {
    Hmat <- build_Hmat_NoError(data_aug)
    
  } else if (error_model == "constanterror") {
    sd_xy    <- 10     # horizontal SD
    sd_depth <- 25     # depth SD
    Hmat <- build_Hmat_ConstantError(data_aug, sd_xy, sd_depth)
    
  } else if (error_model == "linearerror") {
    var0_xy  <- 100    # baseline variance
    var1_xy  <- 2      # HDOP scaling
    sd_depth <- 25
    Hmat <- build_Hmat_LinearError(data_aug, var0_xy, var1_xy, sd_depth)
  }
  
  #### INITIAL STATE
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
  
  #### FILTER
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
# CHOOSE ERROR MODEL HERE
############################################################

error_model_choice <- "linearerror"   # or "noerror" or "constanterror"

############################################################
# FIT MODEL
############################################################

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood(p, aug, error_model = error_model_choice),
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
# BUILD H MATRIX (MATCH error_model_choice)
############################################################

if (error_model_choice == "noerror") {
  
  Hmat <- build_Hmat_NoError(aug)
  
} else if (error_model_choice == "constanterror") {
  
  sd_xy    <- 10
  sd_depth <- 25
  Hmat <- build_Hmat_ConstantError(aug, sd_xy, sd_depth)
  
} else if (error_model_choice == "linearerror") {
  
  var0_xy  <- 100
  var1_xy  <- 2
  sd_depth <- 25
  Hmat <- build_Hmat_LinearError(aug, var0_xy, var1_xy, sd_depth)
}

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

## XY, Depth vs Time, X vs Time, Y vs Time, 3D plot

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










