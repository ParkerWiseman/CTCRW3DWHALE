##########################
## Make Kalman matrices ##
##########################
makeT <- function(beta, dt) {
    ebt <- exp(-beta*dt)
    
    Tmat <- matrix(0, 4, 4)
    Tmat[1,1] <- 1
    Tmat[3,3] <- 1
    Tmat[1,2] <- (1 - ebt) / beta
    Tmat[3,4] <- Tmat[1,2]
    Tmat[2,2] <- ebt
    Tmat[4,4] <- Tmat[2,2]
    
    return(Tmat)
}

makeQ <- function(beta, sigma, dt) {
    ebt <- exp(-beta*dt)
    e2bt <- exp(-2*beta*dt)
    
    Qmat <- matrix(0, 4, 4)
    Qmat[1,1] <- (sigma/beta)^2 * (dt + (1-e2bt)/(2*beta) - 2*(1-ebt)/beta)
    Qmat[3,3] <- Qmat[1,1]
    Qmat[1,2] <- sigma^2/(2*beta^2) * (1 - 2*ebt + e2bt)
    Qmat[2,1] <- Qmat[1,2]
    Qmat[3,4] <- Qmat[1,2]
    Qmat[4,3] <- Qmat[1,2]
    Qmat[2,2] <- sigma^2/(2*beta) * (1 - e2bt)
    Qmat[4,4] <- Qmat[2,2]
    
    return(Qmat)
}