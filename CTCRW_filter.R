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