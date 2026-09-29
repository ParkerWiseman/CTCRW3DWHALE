
#' Smoother for CTCRW model
#' 
#' Check out https://github.com/NMML/crawl/blob/master/src/CTCRWPREDICT.cpp
#' Also see Section 4.4.1-4.4.4 of Durbin & Koopman
#' 
#' @param times Vector of numeric times of observations
#' @param obs Matrix with two columns: x and y coordinates
#' @param tau Vector of values of tau parameter (same length as times)
#' @param nu Vector of values of nu parameter (same length as times)
#' @param a0 Initial state mean in Kalman filter
#' @param P0 Initial state covariance in Kalman filter
#' @param H_array Array of measurement error covariance matrices (each layer
#' is a 2-by-2 covariance matrix, and there is one layer for each position
#' in obs)
#' 
#' @return A list of various quantities needed in the smoother and simulation
#' smoother
smoother_ctcrw <- function(times, obs, tau, nu, a0, P0, H_array) {
    obs <- as.matrix(obs)
    n <- nrow(obs)
    
    # Get required quantities from Kalman filter
    filt <- filter_ctcrw(times = times, 
                         obs = obs, 
                         tau = tau, 
                         nu = nu, 
                         a0 = a0, 
                         P0 = P0, 
                         H_array = H_array)
    aest <- filt$aest
    Pest <- filt$Pest
    v <- filt$v
    F <- filt$F
    L <- filt$L
    
    # Initialise quantities needed in the smoother
    r <- c(0, 0, 0, 0)
    N <- matrix(0, 4, 4)
    Z <- matrix(0, 2, 4)
    Z[1,1] <- 1
    Z[2,3] <- 1
    pred <- matrix(0, nrow = n, ncol = 4)
    predVar <- array(0, c(4, 4, n))
    
    # Run backward loop
    for(i in (n+1):2) {
        if(any(is.na(obs[i-1,])) | F[1,1,i-1] * F[2,2,i-1] == 0) {
            r <- t(L[,,i-1]) %*% r
            N <- t(L[,,i-1]) %*% N %*% L[,,i-1]
        } else {
            r <- t(Z) %*% solve(F[,,i-1], v[i-1,]) + t(L[,,i-1]) %*% r
            N <- t(Z) %*% solve(F[,,i-1], Z) + t(L[,,i-1]) %*% N %*% L[,,i-1]
        }
        pred[i-1,] <- aest[i-1,] + Pest[,,i-1] %*% r
        predVar[,,i-1] <- Pest[,,i-1] - Pest[,,i-1] %*% N %*% Pest[,,i-1]
    }
    
    return(list(pred = pred, predVar = predVar))
}