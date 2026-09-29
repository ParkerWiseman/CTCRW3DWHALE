


############################################################
# SPERM WHALE — 2D CRAWL CTCRW MODEL
#
# Purpose:
#   Fit a horizontal 2D CTCRW model using the crawl package
#   and compare its estimates with the custom 3D CTCRW model.
#
# Dataset:
#   sperm_whale_processed.csv
#
# Important:
#   The original dataset contains approximately 4468 rows,
#   but only 731 rows have valid longitude/latitude positions.
#
#   We KEEP ALL 4468 timestamps.
#   Rows without geographic coordinates remain NA in the
#   projected x/y coordinates.
#
#   This allows CRAWL to retain the temporal structure of
#   the complete dataset rather than simply deleting the
#   missing-location observations.
#
# CRAWL estimates:
#   Horizontal CTCRW movement parameters only:
#       beta
#       sigma
#       autocorrelation time = 1 / beta
#       RMS speed
#       mean speed
#
# It does NOT estimate the vertical component.
############################################################


############################################################
# 1. LOAD PACKAGES
############################################################

library(dplyr)
library(lubridate)
library(sf)
library(crawl)
library(ggplot2)


############################################################
# 2. READ DATA
############################################################

whale <- read.csv(
  "sperm_whale_processed.csv",
  stringsAsFactors = FALSE
)

cat("\n")
cat("====================================================\n")
cat("INITIAL DATASET\n")
cat("====================================================\n")

cat(
  "Rows:",
  nrow(whale),
  "\n"
)

cat(
  "Columns:",
  paste(names(whale), collapse = ", "),
  "\n"
)


############################################################
# 3. CHECK REQUIRED COLUMNS
############################################################

required_columns <- c(
  "time",
  "lon",
  "lat",
  "hdop"
)

missing_columns <- setdiff(
  required_columns,
  names(whale)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "The following required columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


############################################################
# 4. CONVERT TIME
############################################################

whale$time <- ymd_hms(
  whale$time,
  tz = "UTC"
)

if (any(is.na(whale$time))) {
  
  stop(
    "Some timestamps could not be converted to POSIXct."
  )
}


############################################################
# 5. SORT BY TIME
############################################################

whale <- whale %>%
  arrange(time)


############################################################
# 6. REMOVE DUPLICATE TIMESTAMPS
############################################################

n_before_duplicates <- nrow(whale)

whale <- whale[
  !duplicated(whale$time),
]

n_removed_duplicates <-
  n_before_duplicates - nrow(whale)

cat(
  "Duplicate timestamps removed:",
  n_removed_duplicates,
  "\n"
)


############################################################
# 7. REPORT MISSING DATA
############################################################

cat("\n")
cat("====================================================\n")
cat("MISSING DATA CHECK\n")
cat("====================================================\n")

cat(
  "Total rows:",
  nrow(whale),
  "\n"
)

cat(
  "Missing longitude:",
  sum(is.na(whale$lon)),
  "\n"
)

cat(
  "Missing latitude:",
  sum(is.na(whale$lat)),
  "\n"
)

cat(
  "Missing HDOP:",
  sum(is.na(whale$hdop)),
  "\n"
)


############################################################
# 8. CREATE DATASET WITH VALID GEOGRAPHIC POSITIONS
############################################################
#
# IMPORTANT:
#
# We do NOT delete the rows with missing coordinates from
# the main dataset.
#
# We only use rows with valid lon/lat to calculate projected
# coordinates.
#
############################################################

valid_coords <- whale %>%
  filter(
    !is.na(lon),
    !is.na(lat)
  )

cat("\n")
cat(
  "Valid geographic positions:",
  nrow(valid_coords),
  "\n"
)


############################################################
# 9. CHECK THAT VALID COORDINATES EXIST
############################################################

if (nrow(valid_coords) == 0) {
  
  stop(
    "No valid longitude/latitude positions were found."
  )
}


############################################################
# 10. DETERMINE UTM ZONE
############################################################

mean_lon <- mean(
  valid_coords$lon,
  na.rm = TRUE
)

mean_lat <- mean(
  valid_coords$lat,
  na.rm = TRUE
)

utm_zone <- floor(
  (mean_lon + 180) / 6
) + 1

if (mean_lat >= 0) {
  
  epsg_code <- 32600 + utm_zone
  
} else {
  
  epsg_code <- 32700 + utm_zone
  
}


cat("\n")
cat("====================================================\n")
cat("UTM PROJECTION\n")
cat("====================================================\n")

cat(
  "Mean longitude:",
  mean_lon,
  "\n"
)

cat(
  "Mean latitude:",
  mean_lat,
  "\n"
)

cat(
  "UTM zone:",
  utm_zone,
  "\n"
)

cat(
  "EPSG:",
  epsg_code,
  "\n"
)


############################################################
# 11. CONVERT VALID POSITIONS TO SF
############################################################

points_sf <- st_as_sf(
  valid_coords,
  coords = c(
    "lon",
    "lat"
  ),
  crs = 4326,
  remove = FALSE
)


############################################################
# 12. PROJECT TO UTM
############################################################

points_sf_utm <- st_transform(
  points_sf,
  crs = epsg_code
)


############################################################
# 13. EXTRACT PROJECTED COORDINATES
############################################################

coords <- st_coordinates(
  points_sf_utm
)

valid_coords$crawl_x <- coords[, 1]

valid_coords$crawl_y <- coords[, 2]


############################################################
# 14. CREATE FULL-LENGTH X/Y COLUMNS
############################################################
#
# These columns have one value for every original
# observation.
#
# Rows without geographic coordinates remain NA.
#
############################################################

whale$crawl_x <- NA_real_

whale$crawl_y <- NA_real_


############################################################
# 15. MATCH PROJECTED POSITIONS BACK TO FULL DATASET
############################################################

match_index <- match(
  valid_coords$time,
  whale$time
)

whale$crawl_x[match_index] <-
  valid_coords$crawl_x

whale$crawl_y[match_index] <-
  valid_coords$crawl_y


############################################################
# 16. VERIFY PROJECTED COORDINATES
############################################################

cat("\n")
cat("====================================================\n")
cat("PROJECTED COORDINATE CHECK\n")
cat("====================================================\n")

cat(
  "Total observations:",
  nrow(whale),
  "\n"
)

cat(
  "Projected X values:",
  sum(!is.na(whale$crawl_x)),
  "\n"
)

cat(
  "Projected Y values:",
  sum(!is.na(whale$crawl_y)),
  "\n"
)

cat(
  "Missing X/Y observations:",
  sum(
    is.na(whale$crawl_x) |
      is.na(whale$crawl_y)
  ),
  "\n"
)


############################################################
# 17. HANDLE HDOP
############################################################
#
# Missing HDOP values are assigned 5.
#
# This follows the previous CRAWL setup.
#
############################################################

whale$hdop[
  is.na(whale$hdop)
] <- 5


############################################################
# 18. CREATE LOCATION ERROR VARIABLES
############################################################

whale$sd_xy <-
  10 * whale$hdop

whale$ln.sd.x <-
  log(whale$sd_xy)

whale$ln.sd.y <-
  log(whale$sd_xy)

whale$error.corr <- 0


############################################################
# 19. CREATE NUMERIC TIME
############################################################

t0 <- min(
  whale$time,
  na.rm = TRUE
)

whale$TimeNum <- as.numeric(
  difftime(
    whale$time,
    t0,
    units = "secs"
  )
)


############################################################
# 20. CREATE DATASET FOR CRAWL
############################################################
#
# Keep ALL observations.
#
# x and y are NA where no geographic position exists.
#
############################################################

crawl_data <- whale %>%
  select(
    x = crawl_x,
    y = crawl_y,
    TimeNum,
    ln.sd.x,
    ln.sd.y,
    error.corr
  )


############################################################
# 21. DISPLAY CRAWL DATA CHECK
############################################################

cat("\n")
cat("====================================================\n")
cat("CRAWL DATA CHECK\n")
cat("====================================================\n")

cat(
  "Rows supplied to CRAWL:",
  nrow(crawl_data),
  "\n"
)

cat(
  "Rows with geographic positions:",
  sum(
    !is.na(crawl_data$x) &
      !is.na(crawl_data$y)
  ),
  "\n"
)

cat(
  "Rows without geographic positions:",
  sum(
    is.na(crawl_data$x) |
      is.na(crawl_data$y)
  ),
  "\n"
)

cat(
  "Time range (seconds):",
  min(crawl_data$TimeNum),
  "to",
  max(crawl_data$TimeNum),
  "\n"
)


############################################################
# 22. DEFINE CRAWL ERROR MODEL
############################################################

err.model <- list(
  
  x = ~ ln.sd.x - 1,
  
  y = ~ ln.sd.y - 1,
  
  rho = ~ error.corr
  
)


############################################################
# 23. INITIAL MOVEMENT PARAMETERS
############################################################
#
# Initial mean speed:
#       1 m/s
#
# Initial autocorrelation time:
#       8 hours
#
############################################################

mean_speed_init <- 1.0

time_scale_init <- 8 * 3600

beta_init <-
  1 / time_scale_init

sigma_init <-
  mean_speed_init *
  sqrt(
    2 * beta_init
  )

theta_init <- c(
  
  log(sigma_init),
  
  log(beta_init)
  
)


############################################################
# 24. FIX LOCATION ERROR PARAMETERS
############################################################
#
# The first two parameters correspond to the x/y
# measurement-error terms.
#
# We fix these and estimate sigma and beta.
#
############################################################

fixPar <- c(
  1,
  1,
  NA,
  NA
)


############################################################
# 25. FIT CRAWL MODEL
############################################################

cat("\n")
cat("====================================================\n")
cat("FITTING CRAWL MODEL\n")
cat("====================================================\n")

fit <- crwMLE(
  
  data = crawl_data,
  
  coord = c(
    "x",
    "y"
  ),
  
  Time.name = "TimeNum",
  
  err.model = err.model,
  
  theta = theta_init,
  
  fixPar = fixPar,
  
  attempts = 10
  
)


############################################################
# 26. EXTRACT CRAWL MOVEMENT PARAMETERS
############################################################

cat("\n")
cat("====================================================\n")
cat("CRAWL FITTED PARAMETERS\n")
cat("====================================================\n")


############################################################
# SHOW STRUCTURE OF FITTED PARAMETER VECTOR
############################################################

cat("\n")
cat("Raw fit$par:\n")

print(fit$par)


cat("\n")
cat("Length of fit$par:\n")

print(length(fit$par))


cat("\n")
cat("Names of fit$par:\n")

print(names(fit$par))


############################################################
# CRAWL PARAMETER ORDER
############################################################
#
# For this model:
#
#   1 = ln tau.x 1
#   2 = ln tau.y 1
#   3 = ln sigma 1
#   4 = ln beta 1
#
# The first two parameters are fixed using:
#
#   fixPar = c(1, 1, NA, NA)
#
# Therefore:
#
#   fit$par[3] = ln sigma
#   fit$par[4] = ln beta
#
############################################################


ln_sigma_mle <- fit$par[3]

ln_beta_mle <- fit$par[4]


############################################################
# CONVERT LOG PARAMETERS TO ORIGINAL SCALE
############################################################

sigma_mle <- exp(
  ln_sigma_mle
)

beta_mle <- exp(
  ln_beta_mle
)


############################################################
# AUTOCORRELATION TIME
############################################################

tau_mle <- 1 / beta_mle


############################################################
# RMS SPEED
############################################################

rms_speed_mle <-
  sigma_mle /
  sqrt(
    2 * beta_mle
  )


############################################################
# MEAN SPEED
############################################################

mean_speed_mle <-
  sqrt(
    pi / beta_mle
  ) *
  sigma_mle / 2


############################################################
# PRINT MOVEMENT PARAMETERS
############################################################

cat("\n")
cat("----------------------------------------------------\n")
cat("CRAWL MOVEMENT ESTIMATES\n")
cat("----------------------------------------------------\n")

cat(
  "CRAWL ln(sigma):",
  ln_sigma_mle,
  "\n"
)

cat(
  "CRAWL ln(beta):",
  ln_beta_mle,
  "\n"
)

cat(
  "CRAWL sigma:",
  sigma_mle,
  "\n"
)

cat(
  "CRAWL beta:",
  beta_mle,
  "\n"
)

cat(
  "CRAWL autocorrelation time (seconds):",
  tau_mle,
  "\n"
)

cat(
  "CRAWL autocorrelation time (hours):",
  tau_mle / 3600,
  "\n"
)

cat(
  "CRAWL RMS speed (m/s):",
  rms_speed_mle,
  "\n"
)

cat(
  "CRAWL mean speed (m/s):",
  mean_speed_mle,
  "\n"
)


############################################################
# CHECK CONVERGENCE
############################################################

cat("\n")
cat("----------------------------------------------------\n")
cat("CRAWL OPTIMIZATION INFORMATION\n")
cat("----------------------------------------------------\n")

cat(
  "CRAWL convergence code:",
  fit$convergence,
  "\n"
)

if (!is.null(fit$message)) {
  
  cat(
    "CRAWL optimizer message:",
    fit$message,
    "\n"
  )
  
}


############################################################
# CREATE A READABLE PARAMETER TABLE
############################################################

crawl_parameter_table <- data.frame(
  
  parameter = c(
    "ln_tau_x",
    "ln_tau_y",
    "ln_sigma",
    "ln_beta"
  ),
  
  estimate = fit$par,
  
  fixed = c(
    TRUE,
    TRUE,
    FALSE,
    FALSE
  )
  
)

cat("\n")
cat("----------------------------------------------------\n")
cat("CRAWL PARAMETER TABLE\n")
cat("----------------------------------------------------\n")

print(
  crawl_parameter_table
)




############################################################
# 27. PREDICT AT ORIGINAL OBSERVATION TIMES
############################################################

cat("\n")
cat("====================================================\n")
cat("CRAWL PREDICTIONS — ORIGINAL TIMES\n")
cat("====================================================\n")

pred_original <- crwPredict(
  
  object.crwFit = fit,
  
  predTime = crawl_data$TimeNum
  
)


############################################################
# 28. INSPECT PREDICTION OBJECT
############################################################

cat("\n")
cat("CRAWL prediction columns:\n")

print(
  names(pred_original)
)


############################################################
# 29. IDENTIFY PREDICTED X/Y COLUMNS
############################################################

pred_names <- names(
  pred_original
)

x_candidates <- pred_names[
  grepl(
    "x",
    pred_names,
    ignore.case = TRUE
  )
]

y_candidates <- pred_names[
  grepl(
    "y",
    pred_names,
    ignore.case = TRUE
  )
]

cat("\n")
cat("Possible X prediction columns:\n")

print(
  x_candidates
)

cat("\n")
cat("Possible Y prediction columns:\n")

print(
  y_candidates
)


############################################################
# 30. IDENTIFY SMOOTHED POSITION COLUMNS
############################################################
#
# CRAWL prediction objects can contain several x/y-related
# quantities.
#
# We select the first matching prediction column containing
# "x" and "y" after inspecting the names.
#
############################################################

x_col <- x_candidates[
  grepl(
    "mu|pred|x",
    x_candidates,
    ignore.case = TRUE
  )
][1]

y_col <- y_candidates[
  grepl(
    "mu|pred|y",
    y_candidates,
    ignore.case = TRUE
  )
][1]


############################################################
# 31. CREATE CRAWL PLOT DATA
############################################################

crawl_plot <- data.frame(
  
  TimeNum = pred_original$TimeNum,
  
  smooth_x =
    pred_original[[x_col]],
  
  smooth_y =
    pred_original[[y_col]]
  
)


############################################################
# 32. CHECK PREDICTION DATA
############################################################

cat("\n")
cat("====================================================\n")
cat("PREDICTION CHECK\n")
cat("====================================================\n")

cat(
  "Prediction rows:",
  nrow(crawl_plot),
  "\n"
)

cat(
  "Missing predicted X:",
  sum(is.na(crawl_plot$smooth_x)),
  "\n"
)

cat(
  "Missing predicted Y:",
  sum(is.na(crawl_plot$smooth_y)),
  "\n"
)


############################################################
# 33. SAVE ORIGINAL-TIME PREDICTIONS
############################################################

write.csv(
  
  crawl_plot,
  
  "sperm_whale_CRAWL_smoothed_track_observation_times.csv",
  
  row.names = FALSE
  
)


############################################################
# 34. CREATE REGULAR 30-MINUTE PREDICTION GRID
############################################################

pred_times <- seq(
  
  min(
    crawl_data$TimeNum
  ),
  
  max(
    crawl_data$TimeNum
  ),
  
  by = 1800
  
)


############################################################
# 35. PREDICT CRAWL TRACK EVERY 30 MINUTES
############################################################

cat("\n")
cat("====================================================\n")
cat("CRAWL PREDICTIONS — 30 MINUTE GRID\n")
cat("====================================================\n")

pred_30min <- crwPredict(
  
  object.crwFit = fit,
  
  predTime = pred_times
  
)


############################################################
# 36. IDENTIFY 30-MINUTE X/Y COLUMNS
############################################################

pred30_names <- names(
  pred_30min
)

x30_candidates <- pred30_names[
  grepl(
    "x",
    pred30_names,
    ignore.case = TRUE
  )
]

y30_candidates <- pred30_names[
  grepl(
    "y",
    pred30_names,
    ignore.case = TRUE
  )
]


############################################################
# 37. CREATE 30-MINUTE PLOT DATA
############################################################

x30_col <- x30_candidates[
  grepl(
    "mu|pred|x",
    x30_candidates,
    ignore.case = TRUE
  )
][1]

y30_col <- y30_candidates[
  grepl(
    "mu|pred|y",
    y30_candidates,
    ignore.case = TRUE
  )
][1]


crawl_30min <- data.frame(
  
  TimeNum =
    pred_30min$TimeNum,
  
  smooth_x =
    pred_30min[[x30_col]],
  
  smooth_y =
    pred_30min[[y30_col]]
  
)


############################################################
# 38. SAVE 30-MINUTE PREDICTIONS
############################################################

write.csv(
  
  crawl_30min,
  
  "sperm_whale_CRAWL_smoothed_track_30min.csv",
  
  row.names = FALSE
  
)


############################################################
# 39. PLOT X THROUGH TIME
############################################################

p_x <- ggplot() +
  
  geom_point(
    
    data = crawl_data,
    
    aes(
      x = TimeNum,
      y = x
    ),
    
    alpha = 0.4
    
  ) +
  
  geom_line(
    
    data = crawl_plot,
    
    aes(
      x = TimeNum,
      y = smooth_x
    ),
    
    linewidth = 1
    
  ) +
  
  labs(
    
    title =
      "Sperm Whale — CRAWL X Position",
    
    subtitle =
      "Observed positions and smoothed horizontal trajectory",
    
    x =
      "Time (seconds)",
    
    y =
      "X (m)"
    
  ) +
  
  theme_minimal()


print(p_x)


############################################################
# 40. PLOT Y THROUGH TIME
############################################################

p_y <- ggplot() +
  
  geom_point(
    
    data = crawl_data,
    
    aes(
      x = TimeNum,
      y = y
    ),
    
    alpha = 0.4
    
  ) +
  
  geom_line(
    
    data = crawl_plot,
    
    aes(
      x = TimeNum,
      y = smooth_y
    ),
    
    linewidth = 1
    
  ) +
  
  labs(
    
    title =
      "Sperm Whale — CRAWL Y Position",
    
    subtitle =
      "Observed positions and smoothed horizontal trajectory",
    
    x =
      "Time (seconds)",
    
    y =
      "Y (m)"
    
  ) +
  
  theme_minimal()


print(p_y)


############################################################
# 41. 2D CRAWL SMOOTHED TRACK
############################################################

p_track <- ggplot() +
  
  geom_point(
    
    data = crawl_data,
    
    aes(
      x = x,
      y = y
    ),
    
    alpha = 0.4
    
  ) +
  
  geom_path(
    
    data = crawl_plot,
    
    aes(
      x = smooth_x,
      y = smooth_y
    ),
    
    linewidth = 1
    
  ) +
  
  labs(
    
    title =
      "Sperm Whale — CRAWL Smoothed Track",
    
    subtitle =
      "Observed positions and smoothed horizontal trajectory",
    
    x =
      "X (m)",
    
    y =
      "Y (m)"
    
  ) +
  
  theme_minimal()


print(p_track)


############################################################
# 42. SAVE CRAWL MODEL
############################################################

saveRDS(
  
  fit,
  
  "sperm_whale_CRAWL_fit.rds"
  
)


############################################################
# 43. FINAL SUMMARY
############################################################

cat("\n")
cat("====================================================\n")
cat("CRAWL ANALYSIS COMPLETE\n")
cat("====================================================\n")

cat(
  "Total observations:",
  nrow(crawl_data),
  "\n"
)

cat(
  "Observed geographic positions:",
  sum(
    !is.na(crawl_data$x) &
      !is.na(crawl_data$y)
  ),
  "\n"
)

cat(
  "Missing geographic positions:",
  sum(
    is.na(crawl_data$x) |
      is.na(crawl_data$y)
  ),
  "\n"
)

cat(
  "CRAWL beta:",
  beta_mle,
  "\n"
)

cat(
  "CRAWL sigma:",
  sigma_mle,
  "\n"
)

cat(
  "CRAWL autocorrelation time (hours):",
  tau_mle / 3600,
  "\n"
)

cat(
  "CRAWL RMS speed (m/s):",
  rms_speed_mle,
  "\n"
)

cat(
  "CRAWL mean speed (m/s):",
  mean_speed_mle,
  "\n"
)

cat("\n")

cat(
  "Output files:\n"
)

cat(
  "  sperm_whale_CRAWL_smoothed_track_observation_times.csv\n"
)

cat(
  "  sperm_whale_CRAWL_smoothed_track_30min.csv\n"
)

cat(
  "  sperm_whale_CRAWL_fit.rds\n"
)

cat("\n")

cat("DONE.\n")













