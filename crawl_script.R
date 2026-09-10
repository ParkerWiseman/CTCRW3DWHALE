# import all whales data 
files <- list.files(
  path = "~/Documents/Dalhousie/PhD/Simone Panigada/data/Whales/",
  pattern = "\\.csv$",
  full.names = TRUE
)

library(lubridate)
# 2. create a list with the good colonmes
whales_list <- lapply(files, function(file) {
  df <- read.csv(file)
  df <- df[, c(2, 4, 5, 6, 7, 26, 27, 28, 29)]
  # Convert observation times using lubridate
  df$date <- mdy_hms(df$date, tz = "UTC") # transforme data in UTC
  df <- subset(df, radius > 0 & !is.na(lon) & !is.na(lat))  # Remove missing data or data with no measurement error
  df <- subset(df, !is.na(lat) & !is.na(lon)) # Remove lat and lon with na
  df <- df[!duplicated(df$date), ]  # Remove duplicated times
  return(df)
})

# merge all whalws data in a  1 data.frame
whales_combined <- do.call(rbind, whales_list)
head(whales_combined)

####################
## Project to UTM ##
####################
library(sf)
points_sf <- st_as_sf(whales_combined, coords = c("lon", "lat"), crs = 4326)
points_sf_utm <- st_transform(points_sf, crs = st_crs("+proj=utm +zone=31")) # Save the projected sf object

# Extract coordinates from it 
points_sf_utm$x <- st_coordinates(points_sf_utm)[,1]
points_sf_utm$y <- st_coordinates(points_sf_utm)[,2]



#####################
## Fit crawl model for 1 whale ##
#####################
data <- subset(points_sf_utm, id== "232682") # Select only on whale 

library(crawl)
# Add measurement error columns
data <- cbind(data, argosDiag2Cov(Major = data$major, 
                                  Minor = data$minor, 
                                  Orientation = data$orientation)[,1:3])

# Error model (based on Argos ellipse information)
err.model <- list(x = ~ln.sd.x - 1,
                  y = ~ln.sd.y - 1,
                  rho = ~ error.corr)
# Check order of model parameters
displayPar(err.model = err.model, data = data, fixPar = c(1, 1, NA, NA))

# Initial parameters for crawl based on expected average speed (in meters/hour)
# and expected time scale of autocorrelation (in hours)
mean_speed <- 5000 #when i change the initial value, the estimate parameter doesm't change
time_scale <- 48
beta <- 3/time_scale
sigma <- 2*sqrt(beta)*mean_speed/sqrt(pi)
par_crawl <- c(log(sigma), log(beta))

# Fit model
fit <- crwMLE(data = data,
              err.model = err.model,
              fixPar = c(1, 1, NA, NA),
              theta = par_crawl,
              Time.name = "date",
              attempts = 10)

# Estimated movement parameters
beta_mle <- exp(fit$par[4])
sigma_mle <- exp(fit$par[3])
mean_speed_mle <- sqrt(pi/beta_mle) * sigma_mle/2
time_scale_mle <- 3/beta_mle

####################################
## Predict regular smoothed track ##
####################################
# Time grid for the regularisation (here chosen to be hourly)
pred_times <- seq(data$date[1], data$date[nrow(data)], by = "30 min")

# Predict regular track
pred <- crwPredict(object.crwFit = fit, predTime = pred_times)

