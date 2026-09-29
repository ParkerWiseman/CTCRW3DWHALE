
makeZ_vec <- function(is_na, ndim) {
    n <- length(is_na)
    Z <- array(0, c(ndim, 2*ndim, n))
    for(i in 1:ndim) {
        Z[i, 2*(i-1)+1, !is_na] <- 1
    }
    return(Z)
}
makeT_vec <- function(beta, dt, ndim) {
    n <- length(beta)
    T <- array(0, c(2*ndim, 2*ndim, n))
    for(i in 1:ndim) {
        T[2*(i-1)+1, 2*(i-1)+1,] <- 1
        T[2*(i-1)+1, 2*i,] <- (1 - exp(-beta * dt)) / beta
        T[2*i, 2*i,] <- exp(-beta * dt)
    }
    return(T)
}
makeQ_vec <- function(beta, sigma, dt, ndim) {
    n <- length(beta)
    Q <- array(0, c(2*ndim, 2*ndim, n))
    for(i in 1:ndim) {
        Q[2*(i-1)+1, 2*(i-1)+1,] <- (sigma/beta)^2 * 
            (dt - 2/beta*(1-exp(-beta*dt)) + 1/(2*beta)*(1-exp(-2*beta*dt)))
        Q[2*(i-1)+1, 2*i,] <- sigma^2/(2*beta*beta) * 
            (1 - 2*exp(-beta*dt) + exp(-2*beta*dt))
        Q[2*i, 2*(i-1)+1,] <- Q[2*(i-1)+1, 2*i,]
        Q[2*i, 2*i,] <- sigma^2/(2*beta) * (1-exp(-2*beta*dt))
    }
    return(Q)
}
makeB_vec <- function(beta, dt, ndim) {
    n <- length(beta)
    B <- array(0, c(2*ndim, ndim, n))
    for(i in 1:ndim) {
        B[2*(i-1)+1, i,] <- dt - (1 - exp(-beta*dt))/beta
        B[2*i, i,] <- 1 - exp(-beta*dt)
    }
    return(B)
}

#' Simulation smoother for CTCRW model
#' 
#' See Section 4.9 of Durbin & Koopman
simulation_smoother_ctcrw <- function(times, obs, tau, nu, a0, P0, H_array) {
    y <- as.matrix(obs)
    n <- nrow(y)
    
    # Unpack SDE parameters
    mu <- matrix(0, n, 2) # no drift for now
    beta <- 1/tau
    sigma <- 2 * nu / sqrt(pi/beta)
    dt <- c(diff(times), NA)
    ndim <- ncol(mu)
    
    # Initialise simulated quantities
    alpha_plus <- matrix(0, nrow = n+1, ncol = 2*ndim) # simulated states
    alpha_plus[1,] <- mgcv::rmvn(n = 1, mu = a0, V = P0)
    y_plus <- matrix(0, nrow = n, ncol = ndim) # simulated observations
    # Initialise filtered states (for observations and simulations)
    a <- matrix(0, nrow = n+1, ncol = 2*ndim)
    a[1,] <- a0
    a_plus <- matrix(0, nrow = n+1, ncol = 2*ndim)
    a_plus[1,] <- a0
    # Initialise state covariance
    P <- array(0, c(2*ndim, 2*ndim, n+1))
    P[,,1] <- P0
    
    # Quantities that need to be saved for the backward loop
    v <- matrix(0, n, ndim)
    v_plus <- matrix(0, n, ndim)
    L <- array(0, c(2*ndim, 2*ndim, n))
    F <- array(0, c(ndim, ndim, n))
    K <- array(0, c(2*ndim, ndim, n))
    
    # Kalman matrices
    Z <- makeZ_vec(is_na = (rowSums(is.na(y)) > 0), ndim = ndim)
    T <- makeT_vec(beta = beta, dt = dt, ndim = ndim)
    Q <- makeQ_vec(beta = beta, sigma = sigma, dt = dt, ndim = ndim)
    B <- makeB_vec(beta = beta, dt = dt, ndim = ndim)
    
    # Replace NA by abitrary numeric value because NA*0 = NA rather than
    # 0 in R
    y[is.na(y)] <- 0
    
    ##################
    ## Forward loop ##
    ##################
    for(i in 1:n) {
        H <- H_array[,,i]
        
        # Matrices  that are shared between y and y_plus filters
        F[,,i] <- Z[,,i] %*% P[,,i] %*% t(Z[,,i]) + H
        K[,,i] <- T[,,i] %*% P[,,i] %*% t(Z[,,i]) %*% solve(F[,,i])
        L[,,i] <- T[,,i] - K[,,i] %*% Z[,,i]
        P[,,i+1] <- T[,,i] %*% P[,,i] %*% t(L[,,i]) + Q[,,i]
        
        # Forward simulation (the random part)
        mu_y_plus <- as.vector(Z[,,i] %*% alpha_plus[i,])
        y_plus[i,] <- mgcv::rmvn(n = 1, mu = mu_y_plus, V = H)
        mu_alpha_plus <- as.vector(T[,,i] %*% alpha_plus[i,] + B[,,i] %*% mu[i,])
        alpha_plus[i+1,] <- mgcv::rmvn(n = 1, mu = mu_alpha_plus, V = Q[,,i])
        
        # Kalman filter (observed data)
        v[i,] <- y[i,] - Z[,,i] %*% a[i,]
        a[i+1,] <- T[,,i] %*% a[i,] + K[,,i] %*% v[i,] + 
            B[,,i] %*% mu[i,]
        
        # Kalman filter (simulated data)
        v_plus[i,] <- y_plus[i,] - Z[,,i] %*% a_plus[i,]
        a_plus[i+1,] <- T[,,i] %*% a_plus[i,] + K[,,i] %*% v_plus[i,] + 
            B[,,i] %*% mu[i,]
    }
    
    # Initialise smoothing residuals
    r <- matrix(0, nrow = n+1, ncol = 2*ndim)
    r_plus <- matrix(0, nrow = n+1, ncol = 2*ndim)
    # Initialise smoothed states
    alpha_hat <- matrix(0, n+1, 2*ndim)
    alpha_plus_hat <- matrix(0, n+1, 2*ndim)
    
    # The last layer of L doesn't exist (NA), and fixing it to zero makes
    # it cancel out in the next loop
    L[is.na(L)] <- 0
    
    #####################
    ## Backward smooth ##
    #####################
    for(i in n:1) {
        r[i,] <- t(Z[,,i]) %*% solve(F[,,i]) %*% v[i,] + 
            t(L[,,i]) %*% r[i+1,]
        r_plus[i,] <- t(Z[,,i]) %*% solve(F[,,i]) %*% v_plus[i,] + 
            t(L[,,i]) %*% r_plus[i+1,]
        
        alpha_hat[i,] <- a[i,] + P[,,i] %*% r[i,]
        alpha_plus_hat[i,] <- a_plus[i,] + P[,,i] %*% r_plus[i,]
    }
    alpha_tilde <- alpha_plus[1:n,] - alpha_plus_hat[1:n,] + alpha_hat[1:n,]
    
    return(alpha_tilde)
}
