
#' Kalman filter for CTCRW model
#' 
#' Based on https://github.com/NMML/crawl/blob/master/src/CTCRWN2LL.cpp
#' See also Sections 4.3.1-4.3.2 of Durbin & Koopman
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
filter_ctcrw <- function(times, obs, tau, nu, a0, P0, H_array) {
    n <- length(times)
    dt <- diff(times)
    
    # Format as matrix if input is data frame
    obs <- as.matrix(obs)
    
    # Get SDE parameters
    beta <- 1/tau
    sigma <- 2 * nu / sqrt(pi/beta)
    
    # Initialise state and state covariance
    aest <- matrix(0, n, 4)
    aest[1,] <- a0
    Pest <- array(0, c(4, 4, n))
    Pest[,,1] <- P0
    
    # Initialise quantities needed in the Kalman filter    
    Z <- matrix(0, 2, 4)
    Z[1,1] <- 1
    Z[2,3] <- 1
    v <- matrix(0, nrow = n, ncol = 2)
    F <- array(0, c(2, 2, n))
    K <- array(0, c(4, 2, n))
    L <- array(0, c(4, 4, n))
    T_array <- array(0, c(4, 4, n))
    
    # Loop over time
    for(i in 1:(n-1)) {
        # Get Kalman matrices for this iteration
        Tmat <- makeT(beta[i], dt[i])
        T_array[,,i] <- Tmat
        Qmat <- makeQ(beta[i], sigma[i], dt[i])
        
        # Check whether position is missing, and then run Kalman steps
        if(any(is.na(obs[i,]))) {
            aest[i+1,] <- Tmat %*% aest[i,]
            Pest[,,i+1] <- Tmat %*% Pest[,,i] %*% t(Tmat) + Qmat
            L[,,i] <- Tmat
        } else {
            H <- H_array[,,i+1]
            v[i,] <- obs[i,] - Z %*% aest[i,]
            F[,,i] <- Z %*% Pest[,,i] %*% t(Z) + H
            K[,,i] <- Tmat %*% Pest[,,i] %*% t(Z) %*% solve(F[,,i])
            L[,,i] <- Tmat - K[,,i] %*% Z
            aest[i+1,] <- Tmat %*% aest[i,] + K[,,i] %*% v[i,]
            Pest[,,i+1] <- Tmat %*% Pest[,,i] %*% t(L[,,i]) + Qmat
        }
    }
    
    return(list(aest = aest, Pest = Pest, v = v, F = F, K = K, L = L, T_array = T_array))
}