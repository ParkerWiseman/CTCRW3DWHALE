

# Fixes sigmas and only estimates betas

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
# FIXED PROCESS NOISE (FROM SIM / PRIOR)
############################################################

sigma1_fixed <- 30   # horizontal
sigma2_fixed <- 10   # vertical

############################################################
# PARAMETER ESTIMATION: β ONLY, USING UNIFIED LINEAR HDOP MODEL
############################################################

params_start <- c(
  beta1 = log(1),
  beta2 = log(1)
)

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood_whale_beta_only(
    params       = p,
    data_aug     = aug,
    sigma1_fixed = sigma1_fixed,
    sigma2_fixed = sigma2_fixed,
    error_model  = "linearerror"
  ),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)

beta1_hat <- p_hat["beta1"]
beta2_hat <- p_hat["beta2"]

sigma1_hat <- sigma1_fixed
sigma2_hat <- sigma2_fixed

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
# DELTA (MATCHES LIKELIHOOD)
############################################################

delta_raw   <- diff(aug$Time)
delta_fixed <- pmax(delta_raw, 1e-5)
delta       <- c(delta_fixed[1], delta_fixed)

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
# PLOTS
############################################################

ggplot() +
  geom_point(data = aug, aes(x = x, y = y), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = y), color = "red", linewidth = 1) +
  theme_minimal()

ggplot() +
  geom_point(data = aug, aes(x = x, y = depth), color = "blue", alpha = 0.4) +
  geom_path(data = smooth_track, aes(x = x, y = depth), color = "red", linewidth = 1) +
  theme_minimal()

plot_ly() %>%
  add_markers(data = aug, x = ~x, y = ~y, z = ~depth,
              marker = list(color = 'blue', size = 2)) %>%
  add_lines(data = smooth_track, x = ~x, y = ~y, z = ~depth,
            line = list(color = 'red', width = 4)) %>%
  layout(scene = list(
    xaxis = list(title = "X"),
    yaxis = list(title = "Y"),
    zaxis = list(title = "Depth")
  ))






