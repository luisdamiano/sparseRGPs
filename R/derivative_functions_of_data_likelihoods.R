## data likelihood derivatives wrt the latent GP for different likelihoods 
 
## Poisson

## second derivative function
d2log_py_dff_pois <- function(ff, ...)
{
  args <- list(...)
  m <- args$m
  return(-m * exp(ff))
}

## third derivative function
## second derivative function
d3log_py_dff_pois <- function(ff, ...)
{
  args <- list(...)
  m <- args$m
  return(-m * exp(ff))
}

## first derivative function
dlog_py_dff_pois <- function(ff, y, ...)
{
  args <- list(...)
  m <- args$m
  return(-m * exp(ff) + y)
}

## Gradient function of psi = log( p(y|lambda) ) + log( p(ff|theta, x'))
grad_loglik_fn_pois <- function(ff, y, mu, Sigma12, Sigma22, Z, ...)
{
  ## ff are the GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations 
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21
  
  # ff <- as.numeric(ff)
  
  ## (d/dff) p(y|lambda)
  args <- list(...)
  m <- args$m
  d1 <- -m * exp(ff) + y
  
  ## (d/dff) p(ff|theta)
  R <- chol(Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12))
  
  # d2 <- -1/Z * (ff - mu) + ((1/Z) * Sigma12) %*% solve(a = Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12) , 
  #                                                              b = t(Sigma12) %*% (1/Z * (ff - mu)))
  
  d2 <- -1/Z * (ff - mu) + ((1/Z) * Sigma12) %*% solve(a = R, b = solve(a = t(R) , 
                                                       b = t(Sigma12) %*% (1/Z * (ff - mu))))
  
  return(d1 + d2)
}

## Gradient function of psi = log( p(y|lambda) ) + log( p(ff|theta, x'))
grad_loglik_fn_pois_full <- function(ff, y, mu, Sigma, ...)
{
  ## ff are the GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma variance covariance matrix at observed data locations

  ff <- as.numeric(ff)
  
  ## (d/dff) p(y|lambda)
  args <- list(...)
  m <- args$m
  d1 <- -m * exp(ff) + y
  
  ## (d/dff) p(ff|theta)
  d2 <- - solve(a = Sigma, b = (ff - mu))
  
  return(d1 + d2)
}

## Bernoulli

## second derivative function
d2log_py_dff_bern <- function(ff, ...)
{
  args <- list(...)
  y <- args$y
  
  ## compute probability vector pi(ff)
  # pi_ff <- (1 - 1e-4) * my_logistic(x = ff) + (1e-4) * 0.5
  pi_ff <- my_logistic(x = ff)
  
  ## compute first derivative vector of link wrt ff
  dpi_dff <- pi_ff * (1 - pi_ff)
  
  ## compute second derivative vector of link wrt ff
  d2pi_dff2 <- dpi_dff * (1 - 2 * pi_ff) 
  
  return(
    # (y / pi_ff - 1 / (1 - pi_ff) + y / (1 - pi_ff)) * d2pi_dff2 +
    #   (dpi_dff)^2 * (-y / (pi_ff^2) - 1 / ((1 - pi_ff)^2) - y / ((1 - pi_ff)^2) )
    (1 - 2*pi_ff) * (y - pi_ff) - 
      (y*(1 - pi_ff)^2 + pi_ff^2 + y * pi_ff^2)
  )
  
}

## second derivative function
d3log_py_dff_bern <- function(ff, ...)
{
  args <- list(...)
  y <- args$y
  
  ## compute probability vector pi(ff)
  pi_ff <- my_logistic(x = ff)
  # pi_ff <- (1 - 1e-4) * my_logistic(x = ff) + (1e-4) * 0.5
  
  
  ## compute first derivative vector of link wrt ff
  dpi_dff <- pi_ff * (1 - pi_ff)
  
  ## compute second derivative vector of link wrt ff
  d2pi_dff2 <- dpi_dff * (1 - 2 * pi_ff) 
  
  ## compute third derivative vector of link wrt ff
  d3pi_dff3 <- d2pi_dff2 * (1 - 2 * pi_ff) +
    dpi_dff * (-2 * dpi_dff)
  
  return(
    # dpi_dff * d2pi_dff2 * (-y / (pi_ff^2) - 1/((1 - pi_ff)^2) + y/( (1 - pi_ff)^2) ) +
    #   d3pi_dff3 * (y / pi_ff - 1 / (1 - pi_ff) + y / (1 - pi_ff)) +
    #   dpi_dff^3 * (2 * y / (pi_ff^3) - 2 /( (1 - pi_ff)^3 ) - 2 * y / ((1 - pi_ff)^3)) +
    #   2 * dpi_dff * d2pi_dff2 * (-y / (pi_ff^2) - 1 / ((1 - pi_ff)^2) - y / ((1 - pi_ff)^2))
    -2 * dpi_dff * (y - pi_ff) - dpi_dff * (1 - 2 * pi_ff) - 
      (2 * y * (1 - pi_ff) * (-dpi_dff) + 2 * pi_ff * dpi_dff + 2 * y * pi_ff * dpi_dff)
  )
  
}

## testing the second and third derivative functions
# y <- 1
# f <- seq(from = -1, to = 1, by = 1e-4)
# f_plus <- f + 1e-4 / 2
# f_minus <- f - 1e-4 / 2
# 
# d2s <- numeric()
# d3s <- numeric()
# for(i in 1:length(f))
# {
#   d3s[i] <- d3log_py_dff_bern(ff = f[i], y = y)
#   d2s[i] <- (d2log_py_dff_bern(ff = f_plus[i], y = y) - d2log_py_dff_bern(ff = f_minus[i], y = y)) / (1e-4)
# }
# plot(d2s, y = d3s)
# abline(a = 0, b = 1)

## first derivative function
dlog_py_dff_bern <- function(ff, y, ...)
{
  args <- list(...)
  
  ## compute probability vector pi(ff)
  pi_ff <- my_logistic(x = ff)
  # pi_ff <- (1 - 1e-4) * my_logistic(x = ff) + (1e-4) * 0.5
  
  
  ## compute first derivative vector of link wrt ff
  dpi_dff <- pi_ff * (1 - pi_ff)
  
  ## compute second derivative vector of link wrt ff
  d2pi_dff2 <- dpi_dff * (1 - 2 * pi_ff) 
  
  return(
    # (y / pi_ff - 1 / (1 - pi_ff) + y / (1 - pi_ff)) * dpi_dff 
    y * (1 - pi_ff) - pi_ff + y * pi_ff
  )
}

## Gradient function of psi = log( p(y|lambda) ) + log( p(ff|theta, x'))
grad_loglik_fn_bern <- function(ff, y, mu, Sigma12, Sigma22, Z, ...)
{
  ## ff are the GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations 
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21
  
  # ff <- as.numeric(ff)
  ## compute probability vector pi(ff)
  # pi_ff <- (1 - 1e-4) * my_logistic(x = ff) + (1e-4) * 0.5
  pi_ff <- my_logistic(x = ff)
  
  ## compute first derivative vector of link wrt ff
  dpi_dff <- pi_ff * (1 - pi_ff)
  
  ## compute second derivative vector of link wrt ff
  d2pi_dff2 <- dpi_dff * (1 - 2 * pi_ff) 
  
  ## (d/dff) p(y|lambda)
  args <- list(...)
  # d1 <- (y / pi_ff - 1 / (1 - pi_ff) + y / (1 - pi_ff)) * dpi_dff
  d1 <- y * (1 - pi_ff) - pi_ff + y * pi_ff
  
  ## (d/dff) p(ff|theta)
  R <- chol(Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12))
  
  # d2 <- -1/Z * (ff - mu) + ((1/Z) * Sigma12) %*% solve(a = Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12) , 
  #                                                              b = t(Sigma12) %*% (1/Z * (ff - mu)))
  
  d2 <- -1/Z * (ff - mu) + ((1/Z) * Sigma12) %*% solve(a = R, b = solve(a = t(R) , 
                                                                        b = t(Sigma12) %*% (1/Z * (ff - mu))))
  
  if(any(is.nan(d1 + d2)))
  {
    # print("d1 is...")
    # print(d1)
    # 
    # print("d2 is...")
    # print(d2)
    print("min and max ff is...")
    print(c(min(ff), max(ff)))
  }
  
  return(d1 + d2)
}

## Gradient function of psi = log( p(y|lambda) ) + log( p(ff|theta, x'))
grad_loglik_fn_bern_full <- function(ff, y, mu, Sigma, ...)
{
  ## ff are the GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma variance covariance matrix at observed data locations
  
  ff <- as.numeric(ff)
  
  ## compute probability vector pi(ff)
  pi_ff <- my_logistic(x = ff)
  
  ## compute first derivative vector of link wrt ff
  dpi_dff <- pi_ff * (1 - pi_ff)
  
  ## compute second derivative vector of link wrt ff
  d2pi_dff2 <- dpi_dff * (1 - 2 * pi_ff) 
  
  ## (d/dff) p(y|lambda)
  args <- list(...)
  # d1 <- (y / pi_ff - 1 / (1 - pi_ff) + y / (1 - pi_ff)) * dpi_dff
  d1 <- y * (1 - pi_ff) - pi_ff + y * pi_ff
  
  ## (d/dff) p(ff|theta)
  d2 <- - solve(a = Sigma, b = (ff - mu))
  
  return(d1 + d2)
}