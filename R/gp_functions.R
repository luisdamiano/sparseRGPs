## function to get the conditional mean of a normal distribution 
normal_cond_mean <- function(y, x, x_pred, mu, sigma)
{
  ## mu is the same length as (x, x_pred)
  ## sigma is the variance covariance matrix with the rows and columns 
    ##  ordered corresponding to (x,x_pred)
  temp <- nrow(x) + nrow(x_pred)
  
  sigma21 <- sigma[(nrow(x) + 1):temp,1:nrow(x)]
  sigma11 <- sigma[1:nrow(x), 1:nrow(x)]
  return(mu[(nrow(x) + 1):temp] + sigma21 %*% qr.solve(sigma11) %*% (y - mu[1:nrow(x)]))
}

## function to get the conditional variances of a normal distribution 
normal_cond_var <- function(x, x_pred, sigma)
{
  ## sigma is the variance covariance matrix with the rows and columns 
    ##  ordered corresponding to (x,x_pred)
  temp <- nrow(x) + nrow(x_pred)
  
  sigma21 <- sigma[(nrow(x) + 1):temp,1:nrow(x)]
  sigma11 <- sigma[1:nrow(x), 1:nrow(x)]
  sigma22 <- sigma[(nrow(x) + 1):temp,(nrow(x) + 1):temp]
  sigma12 <- sigma[1:nrow(x), (nrow(x) + 1):temp]
  return(sigma22 - sigma21 %*% solve(sigma11) %*% sigma12)
}

## function to scale a matrix based on the column means and standard deviations
my_scale <- function(x, mean_vec, sd_vec)
{
  mean_mat <- matrix(rep(x = mean_vec, times = nrow(x)), 
                     ncol = ncol(x), nrow = nrow(x), byrow = TRUE)
  sd_mat <- matrix(rep(x = sd_vec, times = nrow(x)), 
                   ncol = ncol(x), nrow = nrow(x), byrow = TRUE)
  
  return((x - mean_mat) / sd_mat)
}

# obj <- function(par, y,xy, u, xu, cov_fun, tau)
# {
#   ## remove y values where the x values for y and u overlap
#   sigma <- exp(par[1])
#   l <- exp(par[2])
#   cov_par <- list("sigma" = sigma, "l" = l)
#   
#   temp <- cbind(y,xy)
#   temp <- temp[!(xy %in% xu),]
#   y <- temp[,1]
#   xy <- temp[,2]
# 
#   ## compute the covariance matrix of the active points
#   Kuu <- make_cov_mat(x = xu, x_pred = numeric(), cov_fun = cov_fun, cov_par = cov_par, tau = tau)
# 
#   ## compute Kyu
#   Kyu <- matrix(nrow = length(y), ncol = length(u))
#   for(i in 1:length(y))
#   {
#     for(j in 1:length(u))
#     {
#       Kyu[i,j] <- cov_fun(x1 = xy[i], x2 = xu[j], cov_par)
#     }
#   }
# 
#   ## compute conditional mean and variance for the y values
#   Kuu_inv <- qr.solve(a = Kuu)
#   cond_mean <- Kyu %*% Kuu_inv %*% u
#   cond_var <- sigma^2 + tau^2 - diag(Kyu %*% Kuu_inv %*% t(Kyu))
# 
#   ## compute the joint log pdf values
#   lpdf_y_given_u <- sum(dnorm(x = y, mean = as.numeric(cond_mean), sd = sqrt(cond_var), log = TRUE))
#   lpdf_u <- mvtnorm::dmvnorm(x = u, sigma = Kuu, log = TRUE)
#   lpdf <- lpdf_y_given_u + lpdf_u
#   return(-lpdf)
# }

## function to estimate hyperparameter values via maximum likelihood
# mlgp_low_rank <- function(y, xy, u, xu, sigma, l, tau = 1e-6, cov_fun, obj_fun)
# {
# 
#   temp <- optim(fn = obj_fun, par = c(sigma,l), y = y, xy = xy, u = u, xu = xu, cov_fun = cov_fun, tau = tau)
# 
#   return(temp)
# 
# }


## get the conditional distribution and joint pdf values for the low rank GP
# low_rank_cond_dsn <- function(y,xy, u, xu, sigma, l, tau = 1e-5, cov_fun)
# {
#   ## remove y values where the x values for y and u overlap
#   # temp <- cbind(y,xy)
#   # temp <- temp[!(xy %in% xu),]
#   # y <- temp[,1]
#   # xy <- temp[,2]
# 
#   cov_par <- list("sigma" = sigma, "l" = l)
#   
#   ## compute the covariance matrix of the active points
#   Kuu <- make_cov_mat(x = xu, x_pred = numeric(), cov_fun = cov_fun_sqrd_exp, cov_par = cov_par, tau = tau)
# 
#   ## compute Kyu
#   Kyu <- matrix(nrow = length(y), ncol = length(u))
#   for(i in 1:length(y))
#   {
#     for(j in 1:length(u))
#     {
#       Kyu[i,j] <- cov_fun(x1 = xy[i], x2 = xu[j], cov_par)
#     }
#   }
# 
#   ## compute conditional mean and variance for the y values
#   Kuu_inv <- qr.solve(a = Kuu)
#   cond_mean <- Kyu %*% Kuu_inv %*% u
#   cond_var <- sigma^2 + tau^2 - diag(Kyu %*% Kuu_inv %*% t(Kyu))
# 
#   ## compute the joint log pdf values
#   lpdf_y_given_u <- sum(dnorm(x = y, mean = as.numeric(cond_mean), sd = sqrt(cond_var), log = TRUE))
#   lpdf_u <- mvtnorm::dmvnorm(x = u, sigma = Kuu, log = TRUE)
#   lpdf <- lpdf_y_given_u + lpdf_u
# 
#   ## compute 95% intervals for Y|u
#   lower <- as.numeric(cond_mean) - 1.96 * sqrt(cond_var)
#   upper <- as.numeric(cond_mean) + 1.96 * sqrt(cond_var)
# 
#   return(list("est" = as.numeric(cond_mean), "xy" = xy, "lower" = lower, "upper" = upper, "lpdf" = lpdf, "cond_var" = cond_var))
# 
# }


## log pdf of an AR(1) process/ornstein-uhlenbeck GP with EVENLY SPACED OBSERVATIONS
# lognormal_ou_pdf <- function(x, mu, sigma, l)
# {
#   n <- length(x)
#   rho <- exp(-1/l)
# 
#   return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
#          - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
#                                               - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
# }