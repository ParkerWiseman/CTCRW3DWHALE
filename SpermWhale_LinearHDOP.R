


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






