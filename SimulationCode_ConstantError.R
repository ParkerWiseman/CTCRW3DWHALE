
#Constant Error Simulation Experiment.

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

###############################################
# MATRICES
###############################################


###############################################
# TRUE PARAMETERS
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

s1 <- sigma1_true^2
s2 <- sigma1_true^2
s3 <- sigma2_true^2

###############################################
# SIMULATION
###############################################

N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)

Tmat <- makeT_R(beta1_true, beta2_true, dt)
Qmat <- makeQ_R(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

sd_xy    <- 5
sd_depth <- 5

obs_x     <- X[,1] + rnorm(N, 0, sd_xy)
obs_y     <- X[,3] + rnorm(N, 0, sd_xy)
obs_depth <- X[,5] + rnorm(N, 0, sd_depth)

sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

###############################################
# CONSTANT H MATRIX
###############################################



###############################################
# FILTER
###############################################



###############################################
# SMOOTHER
###############################################



###############################################
# NEGATIVE LOG-LIKELIHOOD
###############################################



###############################################
# OPTIMIZATION
###############################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_ConstantError,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)







#################################
#################################
#################################
#################################

#(ignore below code)













# Parameter values are re-obtained.

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

# TRUE PARAMETERS
beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

s1 <- sigma1_true^2
s2 <- sigma1_true^2
s3 <- sigma2_true^2

# TIME GRID
N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

# LATENT STATE: (x, vx, y, vy, depth, vdepth)
X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)



# SIMULATE LATENT PROCESS
Tmat <- makeT_R(beta1_true, beta2_true, dt)
Qmat <- makeQ_R(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

# CONSTANT SMALL MEASUREMENT ERROR
sd_xy     <- 5     # horizontal SD
sd_depth  <- 5     # depth SD

obs_x     <- X[,1] + rnorm(N, 0, sd_xy)
obs_y     <- X[,3] + rnorm(N, 0, sd_xy)
obs_depth <- X[,5] + rnorm(N, 0, sd_depth)

sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])





# KALMAN FILTER
CTCRW_filter <- function(y, Hmat, beta1_vec, beta2_vec,
                         s1, s2, s3, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT_R(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ_R(beta1_vec[i], beta2_vec[i], s1, s2, s3, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    Z_i <- Z_full
    H_i <- diag(Hmat[i,], 3)
    
    v <- y[i,] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,]  <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}

# SMOOTHER
CTCRW_smoother <- function(filter_out, beta1_vec, beta2_vec,
                           s1, s2, s3, delta) {
  
  N <- nrow(filter_out$a_f)
  a_s <- filter_out$a_f
  P_s <- filter_out$P_f
  
  for (i in (N-1):1) {
    
    Tmat <- makeT_R(beta1_vec[i+1], beta2_vec[i+1], delta[i+1])
    
    a_f_i   <- filter_out$a_f[i,]
    a_p_ip1 <- filter_out$a_p[i+1,]
    P_f_i   <- filter_out$P_f[[i]]
    P_p_ip1 <- filter_out$P_p[[i+1]]
    
    J <- tryCatch(P_f_i %*% t(Tmat) %*% solve(P_p_ip1), error=function(e) NULL)
    if (is.null(J)) next
    
    a_s[i,]  <- a_f_i + J %*% (a_s[i+1,] - a_p_ip1)
    P_s[[i]] <- P_f_i + J %*% (P_s[[i+1]] - P_p_ip1) %*% t(J)
  }
  
  list(a_s = a_s, P_s = P_s)
}

# NEGATIVE LOG-LIKELIHOOD
neg_loglikelihood <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  Hmat <- build_Hmat_noerror(nrow(data_aug), sd_xy, sd_depth)
  
  delta <- rep(dt, nrow(data_aug))
  
  s1 <- sigma1^2
  s2 <- sigma1^2
  s3 <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  a <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
  P <- diag(6) * 1e2
  
  filt <- CTCRW_filter(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s1        = s1,
    s2        = s2,
    s3        = s3,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  -filt$ll
}

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
p_hat






##################################################################












# Parameter values are re-obtained.

library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

source("matrices.R")

# TRUE PARAMETERS
beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

s1 <- sigma1_true^2
s2 <- sigma1_true^2
s3 <- sigma2_true^2

# TIME GRID
N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

# LATENT STATE: (x, vx, y, vy, depth, vdepth)
X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)



# SIMULATE LATENT PROCESS
Tmat <- makeT(beta1_true, beta2_true, dt)
Qmat <- makeQ(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

# CONSTANT SMALL MEASUREMENT ERROR
sd_xy     <- 5     # horizontal SD
sd_depth  <- 5     # depth SD

obs_x     <- X[,1] + rnorm(N, 0, sd_xy)
obs_y     <- X[,3] + rnorm(N, 0, sd_xy)
obs_depth <- X[,5] + rnorm(N, 0, sd_depth)

sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])





# KALMAN FILTER
CTCRW_filter <- function(y, Hmat, beta1_vec, beta2_vec,
                         s1, s2, s3, delta, a, P) {
  
  N <- nrow(y)
  
  Z_full <- matrix(0, 3, 6)
  Z_full[1,1] <- 1
  Z_full[2,3] <- 1
  Z_full[3,5] <- 1
  
  a_f <- matrix(NA, N, 6)
  P_f <- vector("list", N)
  a_p <- matrix(NA, N, 6)
  P_p <- vector("list", N)
  
  aest <- a
  Pest <- P
  ll <- 0
  
  for (i in 1:N) {
    
    Tmat <- makeT_R(beta1_vec[i], beta2_vec[i], delta[i])
    Qmat <- makeQ_R(beta1_vec[i], beta2_vec[i], s1, s2, s3, delta[i])
    
    a_pred <- Tmat %*% aest
    P_pred <- Tmat %*% Pest %*% t(Tmat) + Qmat
    
    a_p[i,]  <- as.numeric(a_pred)
    P_p[[i]] <- P_pred
    
    Z_i <- Z_full
    H_i <- diag(Hmat[i,], 3)
    
    v <- y[i,] - Z_i %*% aest
    Fmat <- Z_i %*% Pest %*% t(Z_i) + H_i
    
    invF <- solve(Fmat)
    ll <- ll - 0.5 * (log(det(Fmat)) + t(v) %*% invF %*% v)
    
    K <- Tmat %*% Pest %*% t(Z_i) %*% invF
    
    aest <- a_pred + K %*% v
    Pest <- Tmat %*% Pest %*% t(Tmat - K %*% Z_i) + Qmat
    
    a_f[i,]  <- as.numeric(aest)
    P_f[[i]] <- Pest
  }
  
  list(ll = ll, a_f = a_f, P_f = P_f, a_p = a_p, P_p = P_p)
}

# SMOOTHER
CTCRW_smoother <- function(filter_out, beta1_vec, beta2_vec,
                           s1, s2, s3, delta) {
  
  N <- nrow(filter_out$a_f)
  a_s <- filter_out$a_f
  P_s <- filter_out$P_f
  
  for (i in (N-1):1) {
    
    Tmat <- makeT_R(beta1_vec[i+1], beta2_vec[i+1], delta[i+1])
    
    a_f_i   <- filter_out$a_f[i,]
    a_p_ip1 <- filter_out$a_p[i+1,]
    P_f_i   <- filter_out$P_f[[i]]
    P_p_ip1 <- filter_out$P_p[[i+1]]
    
    J <- tryCatch(P_f_i %*% t(Tmat) %*% solve(P_p_ip1), error=function(e) NULL)
    if (is.null(J)) next
    
    a_s[i,]  <- a_f_i + J %*% (a_s[i+1,] - a_p_ip1)
    P_s[[i]] <- P_f_i + J %*% (P_s[[i+1]] - P_p_ip1) %*% t(J)
  }
  
  list(a_s = a_s, P_s = P_s)
}

# NEGATIVE LOG-LIKELIHOOD
neg_loglikelihood <- function(params, data_aug) {
  
  beta1  <- exp(params["beta1"])
  beta2  <- exp(params["beta2"])
  sigma1 <- exp(params["sigma1"])
  sigma2 <- exp(params["sigma2"])
  
  y <- as.matrix(data_aug[, c("x","y","depth")])
  
  Hmat <- build_Hmat_noerror(nrow(data_aug), sd_xy, sd_depth)
  
  delta <- rep(dt, nrow(data_aug))
  
  s1 <- sigma1^2
  s2 <- sigma1^2
  s3 <- sigma2^2
  
  beta1_vec <- rep(beta1, nrow(data_aug))
  beta2_vec <- rep(beta2, nrow(data_aug))
  
  a <- c(y[1,1], 0, y[1,2], 0, y[1,3], 0)
  P <- diag(6) * 1e2
  
  filt <- CTCRW_filter(
    y         = y,
    Hmat      = Hmat,
    beta1_vec = beta1_vec,
    beta2_vec = beta2_vec,
    s1        = s1,
    s2        = s2,
    s3        = s3,
    delta     = delta,
    a         = a,
    P         = P
  )
  
  -filt$ll
}

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
p_hat




###########################################################




# Parameter values are re-obtained.

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

# TRUE PARAMETERS
beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

s1 <- sigma1_true^2
s2 <- sigma1_true^2
s3 <- sigma2_true^2

# TIME GRID
N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

# LATENT STATE: (x, vx, y, vy, depth, vdepth)
X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)



# SIMULATE LATENT PROCESS
Tmat <- makeT_R(beta1_true, beta2_true, dt)
Qmat <- makeQ_R(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

# CONSTANT SMALL MEASUREMENT ERROR
sd_xy     <- 5     # horizontal SD
sd_depth  <- 5     # depth SD

obs_x     <- X[,1] + rnorm(N, 0, sd_xy)
obs_y     <- X[,3] + rnorm(N, 0, sd_xy)
obs_depth <- X[,5] + rnorm(N, 0, sd_depth)

sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])







params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
p_hat












###########################################










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

###############################################
# MATRICES
###############################################


###############################################
# TRUE PARAMETERS
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

s1 <- sigma1_true^2
s2 <- sigma1_true^2
s3 <- sigma2_true^2

###############################################
# SIMULATION
###############################################

N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)

Tmat <- makeT_R(beta1_true, beta2_true, dt)
Qmat <- makeQ_R(beta1_true, beta2_true, s1, s2, s3, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

sd_xy    <- 5
sd_depth <- 5

obs_x     <- X[,1] + rnorm(N, 0, sd_xy)
obs_y     <- X[,3] + rnorm(N, 0, sd_xy)
obs_depth <- X[,5] + rnorm(N, 0, sd_depth)

sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

###############################################
# CONSTANT H MATRIX
###############################################



###############################################
# FILTER
###############################################



###############################################
# SMOOTHER
###############################################



###############################################
# NEGATIVE LOG-LIKELIHOOD
###############################################



###############################################
# OPTIMIZATION
###############################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_ConstantError,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)









################################
################################











# (using s_horiz and s_vert)


library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
set.seed(123)

source("matrices.R")          # contains makeT() and makeQ()
source("CTCRW_filter.R")      # contains CTCRW_filter1()
source("CTCRW_smoother.R")    # contains CTCRW_smoother1()
source("neg_loglikelihood.R") # we will call the updated 2-sigma version

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30   # horizontal
sigma2_true <- 10   # vertical

s_horiz_true <- sigma1_true^2
s_vert_true  <- sigma2_true^2

###############################################
# SIMULATION (2-sigma model)
###############################################

N  <- 5000
dt <- 15 / (24*60)
time_vec <- seq(0, by = dt, length.out = N)

X <- matrix(0, nrow = N, ncol = 6)
X[1,] <- c(0, 0, 0, 0, -50, 0)

Tmat <- makeT(beta1_true, beta2_true, dt)
Qmat <- makeQ(beta1_true, beta2_true, s_horiz_true, s_vert_true, dt)

for (i in 2:N) {
  X[i,] <- Tmat %*% X[i-1,] + MASS::mvrnorm(1, rep(0,6), Qmat)
}

# Add measurement error
sd_xy    <- 5
sd_depth <- 5

obs_x     <- X[,1] + rnorm(N, 0, sd_xy)
obs_y     <- X[,3] + rnorm(N, 0, sd_xy)
obs_depth <- X[,5] + rnorm(N, 0, sd_depth)

sim_data <- data.frame(
  time  = ymd_hms("2020-01-01 00:00:00") + time_vec*86400,
  x     = obs_x,
  y     = obs_y,
  depth = obs_depth
)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])


###############################################
# OPTIMIZATION
###############################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_ConstantError2,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)






#################################

















# IGNORE ALL ABOVE CODE
# GOOD CODE


# uses s_horiz and s_vert, and simulation function






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
source("simulate_CTCRW_3D.R") # NEW

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

###############################################
# SIMULATION USING THE NEW FUNCTION
###############################################

N  <- 5000
dt <- 15 / (24*60)

sim_data <- simulate_CTCRW_3D(
  N           = N,
  dt          = dt,
  beta1_true  = beta1_true,
  beta2_true  = beta2_true,
  sigma1_true = sigma1_true,
  sigma2_true = sigma2_true
)

# Add constant measurement error
sd_xy    <- 5
sd_depth <- 5

sim_data$x     <- sim_data$x     + rnorm(N, 0, sd_xy)
sim_data$y     <- sim_data$y     + rnorm(N, 0, sd_xy)
sim_data$depth <- sim_data$depth + rnorm(N, 0, sd_depth)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")) + 1,
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

###############################################
# OPTIMIZATION
###############################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = neg_loglikelihood_ConstantError2,
  data_aug = aug,
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)









############################

# August 12, 2026


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

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

###############################################
# SIMULATION USING THE NEW FUNCTION
###############################################

N  <- 5000
dt <- 15 / (24*60)

sim_data <- simulate_CTCRW_3D(
  N           = N,
  dt          = dt,
  beta1_true  = beta1_true,
  beta2_true  = beta2_true,
  sigma1_true = sigma1_true,
  sigma2_true = sigma2_true
)

# Add constant measurement error
sd_xy    <- 5
sd_depth <- 5

sim_data$x     <- sim_data$x     + rnorm(N, 0, sd_xy)
sim_data$y     <- sim_data$y     + rnorm(N, 0, sd_xy)
sim_data$depth <- sim_data$depth + rnorm(N, 0, sd_depth)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

###############################################
# OPTIMIZATION USING UNIFIED LIKELIHOOD
###############################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood(p, aug, error_model = "constanterror"),
  method   = "L-BFGS-B",
  control  = list(trace = 1, maxit = 1000)
)

p_hat <- exp(fit$par)
print(p_hat)


















################################
# Same as above code but with plots




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

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

###############################################
# SIMULATION USING THE NEW FUNCTION
###############################################

N  <- 5000
dt <- 15 / (24*60)

sim_data <- simulate_CTCRW_3D(
  N           = N,
  dt          = dt,
  beta1_true  = beta1_true,
  beta2_true  = beta2_true,
  sigma1_true = sigma1_true,
  sigma2_true = sigma2_true
)

# Add constant measurement error
sd_xy    <- 5
sd_depth <- 5

sim_data$x     <- sim_data$x     + rnorm(N, 0, sd_xy)
sim_data$y     <- sim_data$y     + rnorm(N, 0, sd_xy)
sim_data$depth <- sim_data$depth + rnorm(N, 0, sd_depth)

aug <- sim_data %>%
  mutate(
    Time = as.numeric(difftime(time, min(time), units = "days")),
    orig_index = seq_len(n())
  )

y <- as.matrix(aug[, c("x","y","depth")])

###############################################
# OPTIMIZATION USING UNIFIED CONSTANT-ERROR LIKELIHOOD
###############################################

params_start <- c(
  beta1  = log(1),
  beta2  = log(1),
  sigma1 = log(10),
  sigma2 = log(10)
)

fit <- optim(
  par      = params_start,
  fn       = function(p) neg_loglikelihood(p, aug, error_model = "constanterror"),
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

###############################################
# BUILD CONSTANT-ERROR H MATRIX
###############################################

Hmat <- build_Hmat_ConstantError2(aug, sd_xy, sd_depth)

###############################################
# DELTA (MATCHES LIKELIHOOD)
###############################################

delta_raw   <- diff(aug$Time)
delta_fixed <- pmax(delta_raw, 1e-5)
delta       <- c(delta_fixed[1], delta_fixed)

###############################################
# INITIAL STATE
###############################################

a0 <- c(
  y[1,1], 0,
  y[1,2], 0,
  y[1,3], 0
)
P0 <- diag(6) * 1e2

###############################################
# FILTER
###############################################

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

###############################################
# SMOOTHER
###############################################

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

###############################################
# 3D PLOT (ROBUST VERSION)
###############################################

plot_ly() %>%
  add_trace(
    data = smooth_track,
    x = ~x, y = ~y, z = ~depth,
    type = "scatter3d",
    mode = "lines",
    line = list(color = 'red', width = 6),
    name = "Smoothed"
  ) %>%
  add_trace(
    data = aug,
    x = ~x, y = ~y, z = ~depth,
    type = "scatter3d",
    mode = "markers",
    marker = list(color = 'blue', size = 2),
    name = "Observed"
  ) %>%
  layout(
    title = "3D CTCRW Simulation (Constant Error): Observed (blue) vs Smoothed (red)",
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Depth")
    )
  )
























###################################

# histograms from multiple simulations




library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)

source("matrices.R")
source("CTCRW_filter.R")
source("CTCRW_smoother.R")
source("neg_loglikelihood.R")
source("simulate_CTCRW_3D.R")

set.seed(123)

###############################################
# TRUE PARAMETERS (2-sigma model)
###############################################

beta1_true  <- 2.8
beta2_true  <- 0.8
sigma1_true <- 30
sigma2_true <- 10

###############################################
# MONTE CARLO SETTINGS
###############################################

n_sims <- 15
N      <- 5000
dt     <- 15 / (24*60)

# storage for parameter estimates
estimates <- data.frame(
  beta1  = numeric(n_sims),
  beta2  = numeric(n_sims),
  sigma1 = numeric(n_sims),
  sigma2 = numeric(n_sims)
)

###############################################
# MONTE CARLO LOOP
###############################################

for (s in 1:n_sims) {
  
  # --- Simulate latent CTCRW ---
  sim_data <- simulate_CTCRW_3D(
    N           = N,
    dt          = dt,
    beta1_true  = beta1_true,
    beta2_true  = beta2_true,
    sigma1_true = sigma1_true,
    sigma2_true = sigma2_true
  )
  
  # --- Add constant measurement error ---
  sd_xy    <- 5
  sd_depth <- 5
  
  sim_data$x     <- sim_data$x     + rnorm(N, 0, sd_xy)
  sim_data$y     <- sim_data$y     + rnorm(N, 0, sd_xy)
  sim_data$depth <- sim_data$depth + rnorm(N, 0, sd_depth)
  
  # --- Build augmented dataset ---
  aug <- sim_data %>%
    mutate(
      Time = as.numeric(difftime(time, min(time), units = "days")),
      orig_index = seq_len(n())
    )
  
  # --- Optimization ---
  params_start <- c(
    beta1  = log(1),
    beta2  = log(1),
    sigma1 = log(10),
    sigma2 = log(10)
  )
  
  fit <- optim(
    par      = params_start,
    fn       = function(p) neg_loglikelihood(p, aug, error_model = "constanterror"),
    method   = "L-BFGS-B",
    control  = list(trace = 0, maxit = 1000)
  )
  
  p_hat <- exp(fit$par)
  
  # store estimates
  estimates[s, ] <- p_hat
}

###############################################
# HISTOGRAMS OF PARAMETER ESTIMATES
###############################################

par(mfrow=c(2,2))

hist(estimates$beta1, main="β1 estimates", xlab="beta1", col="skyblue")
abline(v=beta1_true, col="red", lwd=2)

hist(estimates$beta2, main="β2 estimates", xlab="beta2", col="skyblue")
abline(v=beta2_true, col="red", lwd=2)

hist(estimates$sigma1, main="σ1 estimates", xlab="sigma1", col="skyblue")
abline(v=sigma1_true, col="red", lwd=2)

hist(estimates$sigma2, main="σ2 estimates", xlab="sigma2", col="skyblue")
abline(v=sigma2_true, col="red", lwd=2)

par(mfrow=c(1,1))

###############################################
# SUMMARY STATISTICS
###############################################

summary(estimates)

