


make_gamma <- function(
    covariate,
    gamma_intercept,
    gamma_coef,
    df = 6
) {
  
  B <- splines::bs(
    covariate,
    df = df,
    degree = 3
  )
  
  gamma_intercept +
    as.vector(B %*% gamma_coef)
}






