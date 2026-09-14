


# Using crawl code to analyze sperm whale dataset to obtain paramter estimates
# (to see if paramter estimates obtained from my code are reasonable)


############################################################
# CRAWL analysis for sperm whale dataset
############################################################

library(dplyr)
library(lubridate)
library(sf)
library(crawl)

############################################################
# LOAD AND PREPARE DATA
############################################################

whale <- read.csv("sperm_whale_processed.csv")

# Convert timestamps to POSIXct
whale$time <- ymd_hms(whale$time, tz = "UTC")

# Remove missing coordinates
whale <- whale %>% filter(!is.na(lon), !is.na(lat))

# Remove duplicated timestamps
whale <- whale[!duplicated(whale$time), ]

############################################################
# PROJECT TO UTM
############################################################

points_sf <- st_as_sf(whale, coords = c("lon", "lat"), crs = 4326)
points_sf_utm <- st_transform(points_sf, crs = st_crs("+proj=utm +zone=31"))

coords <- st_coordinates(points_sf_utm)
whale$x <- coords[,1]
whale$y <- coords[,2]

############################################################
# SIMPLE MEASUREMENT ERROR MODEL
############################################################

# Replace missing HDOP
whale$hdop[is.na(whale$hdop)] <- 5

# Convert HDOP → SD (meters)
whale$sd_xy <- 10 * whale$hdop

# Log SDs required by CRAWL
whale$ln.sd.x <- log(whale$sd_xy)
whale$ln.sd.y <- log(whale$sd_xy)

# No correlation in measurement error
whale$error.corr <- 0

############################################################
# PREPARE DATA FOR CRAWL
############################################################

# Drop geometry
whale <- whale %>% sf::st_drop_geometry()

# Force POSIXct again (for safety)
whale$time <- as.POSIXct(whale$time, tz = "UTC")

# Convert time to numeric seconds since first fix
t0 <- min(whale$time)
whale$TimeNum <- as.numeric(difftime(whale$time, t0, units = "secs"))

# Keep only columns CRAWL needs
whale <- whale %>% 
  select(x, y, TimeNum, ln.sd.x, ln.sd.y, error.corr)

############################################################
# ERROR MODEL
############################################################

err.model <- list(
  x   = ~ ln.sd.x - 1,
  y   = ~ ln.sd.y - 1,
  rho = ~ error.corr
)

############################################################
# INITIAL VALUES
############################################################

mean_speed <- 1.0      # m/s
time_scale <- 8 * 3600 # seconds (for interpretation)
beta_init  <- 1 / time_scale
sigma_init <- mean_speed * sqrt(2 * beta_init)

theta_init <- c(log(sigma_init), log(beta_init))

############################################################
# FIT CRAWL MODEL
############################################################

fit <- crwMLE(
  data       = whale,
  coord      = c("x", "y"),
  time       = "TimeNum",   # numeric time
  time.scale = 1,           # units: seconds
  err.model  = err.model,
  theta      = theta_init,
  fixPar     = c(1, 1, NA, NA),
  attempts   = 10
)

############################################################
# EXTRACT PARAMETER ESTIMATES
############################################################

beta_mle  <- exp(fit$par["beta"])
sigma_mle <- exp(fit$par["sigma"])

mean_speed_mle <- sqrt(pi / beta_mle) * sigma_mle / 2
time_scale_mle <- 1 / beta_mle   # in seconds

cat("Estimated beta:", beta_mle, "\n")
cat("Estimated sigma:", sigma_mle, "\n")
cat("Estimated mean speed (m/s):", mean_speed_mle, "\n")
cat("Estimated autocorrelation timescale (s):", time_scale_mle, "\n")

############################################################
# SMOOTHED TRACK (REGULAR TIME GRID)
############################################################

pred_times <- seq(min(whale$TimeNum), max(whale$TimeNum), by = 1800)  # 30 min in seconds

pred <- crwPredict(
  object.crwFit = fit,
  predTime      = pred_times
)

head(pred)






