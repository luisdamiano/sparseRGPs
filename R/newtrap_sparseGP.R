## My Newton-Raphson algorithm for finding the mode 
##  of the Laplace approximation to a GP model with non-Gaussian likelihood
#source("gp_functions.R")
# Rcpp::sourceCpp(file = "covariance_functions.cpp")
newtrap_sparseGP <- function(start_vals, 
                             obj_fun, 
                             grad_loglik_fn, 
                             dlog_py_dff, 
                             d2log_py_dff, 
                             maxit = 1000, 
                             tol = 1e-6, 
                             cov_par, 
                             cov_fun, 
                             xy, 
                             xu, 
                             y, 
                             mu, 
                             muu, 
                             delta = 1e-6,
                             ...)
{
  ## start_vals is a vector of starting values of the GP corresponding to the observed data locations
  ## maxit is the maximum number of iterations of the algorithm
  ## tol is the absolute tolerance level which dictates 
  ##      how small a difference the algorithm is allowed to make before stopping
  ##      to the objective function as well as how close the gradient is allowed to be to zero
  ## obj_fun is the objective function
  ## d2log_py_dff is a function that computes the vector of second derivatives 
  ##      of log(p(y|lambda)) wrt ff
  ## grad_loglik_fun is a function that computes the gradient of 
  ##      psi = log(p(y|lambda)) + log(p(ff|theta)) wrt ff
  ## cov_par is a list of the covariance function parameters 
  ## cov_fun is a string specifying either "sqexp" or "exp"
  ## xy are the observed data locations
  ## xu are the unobserved knot locations
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## muu is the mean of the GP at each of the knot locations
  ## ... are argument to be passed to the objective function and the gradient of the log-likelihood 
  
  ## create Sigma12
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xu), sep = "")
    Sigma12 <- make_cov_mat_ardC(x = xy, x_pred = xu, cov_par = cov_par,
                                 cov_fun = cov_fun, delta = delta,
                                 lnames = lnames)
    
    ## create Sigma22
    Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(), 
                                 cov_fun = cov_fun, cov_par = cov_par,
                                 delta = delta, lnames = lnames)
  }
  else{
    Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
    
    ## create Sigma22
    Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)
  }
  
  ## create Z
  Z2 <- solve(a = Sigma22, b = t(Sigma12))
  Z3 <- Sigma12 * t(Z2)
  Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
  Z <- cov_par$sigma^2 + cov_par$tau^2 + delta - Z4
  
  ## initialize vector of objective function values
  ff <- start_vals
  obj_fun_vals <-  numeric()
  # if(any(is.nan(ff)))
  # {
  #   print("some f value is NaN at the start of the Newton-Raphson")
  # }
  obj_fun_vals[1] <- obj_fun(ff = ff, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, Z = Z, mu = mu, ...)
  
  ## run NR updates until convergence criteria is met or the maxit is reached 
  iter <- 1
  
  ## update the iteration counter
  iter <- iter + 1
  
  ## update ff
  W <- as.numeric(d2log_py_dff(ff = ff, y = y, ...))
  grad_psi <- grad_loglik_fn(ff = ff, y = y, mu = mu, Sigma12 = Sigma12, Sigma22 = Sigma22, Z = Z, ...)
  ff <- newtrap_sparseGP_update(ff = ff, W = W, Z = Z, 
                                Sigma12 = Sigma12, 
                                Sigma22 = Sigma22, 
                                grad_psi = grad_psi, 
                                dlog_py_dff = dlog_py_dff, y = y, mu = mu, ...)
  
  # if(any(is.nan(ff)))
  # {
  #   print("some f value is NaN after the first update before while loop of the Newton-Raphson")
  # }
  
  ## compute objective function   
  obj_fun_vals[iter] <- obj_fun(ff = ff, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, Z = Z, mu = mu, ...)
  
  while(iter < maxit && (abs(obj_fun_vals[iter] - obj_fun_vals[iter - 1]) > tol || any(abs(grad_psi) > tol)))
  {
    ## update the iteration counter
    iter <- iter + 1
    
    ## update ff
    W <- d2log_py_dff(ff = ff, y = y,  ...)
    grad_psi <- grad_loglik_fn(ff = ff, y = y, mu = mu, Sigma12 = Sigma12, Sigma22 = Sigma22, Z = Z, ...)
    ff <- newtrap_sparseGP_update(ff = ff, W = W, Z = Z, 
                                  Sigma12 = Sigma12, 
                                  Sigma22 = Sigma22, 
                                  grad_psi = grad_psi, 
                                  dlog_py_dff = dlog_py_dff, y = y, mu = mu,...)    
    # if(iter == 63)
    # {
    #   print("63rd iteration values right before crap hits the fan")
    #   print("printing ff...")
    #   global_ff <<- ff
    #   global_cov_par <<- cov_par
    #   print(ff)
    #   print("printing W...")
    #   print(W)
    #   print("printing grad_psi...")
    #   print(grad_psi)
    #   print("sequence of objective function values...")
    #   print(obj_fun_vals)
    # }
    
    # if(any(is.nan(ff)) || any(abs(ff) == Inf || any(abs(ff) > 50)))
    # {
    #   print("some f value is NaN, inf, or large in while loop of the Newton-Raphson")
    #   print("printing ff...")
    #   print(ff)
    #   print("printing W...")
    #   print(W)
    #   print("printing grad_psi...")
    #   print(grad_psi)
    #   print("sequence of objective function values...")
    #   print(obj_fun_vals)
    # }
    
    ## compute objective function   
    # print(obj_fun(ff = ff, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, Z = Z, mu = mu, ...))
    obj_fun_vals[iter] <- obj_fun(ff = ff, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, Z = Z, mu = mu, ...)
    # if(any(is.nan(obj_fun_vals)))
    # {
    #   print(summary(ff))
    #   print(summary(grad_psi))
    #   print(summary(W))
    # }
  }
  
  ## calculate the posterior mean and variances of the latent GP (f or ff) at observed locations
  ##    as well as at the knot locations xu
  
  ## Compute Z^(-1) Sigma12
  ZSig12 <- (1/Z) * Sigma12
  
  ## Compute (W_inv - Z)^-1
  WmZ_inv <- 1/((1/W) - Z)
  
  ## compute TT = t(Sigma12) %*% ((W_inv - Z)^-1) %*% Sigma12
  TT <- t(Sigma12) %*% (WmZ_inv * Sigma12)
  
  ## compute  mu + Sigma21 %*% Sigma11_inv %*% (fhat - mu)
  ##    the posterior mean of u, the GP at the knot location
  R <- chol(Sigma22 + t(Sigma12) %*% ZSig12)
  
  # u_mean <- muu + as.numeric(t(ZSig12) %*% (ff - mu)) - 
  #   as.numeric( t(Sigma12) %*% (ZSig12 %*% 
  #   solve(a = Sigma22 + t(Sigma12) %*% ZSig12, b = t(ZSig12) %*% (ff - mu))) )
  u_mean <- muu + as.numeric(t(ZSig12) %*% (ff - mu)) - 
    as.numeric( t(Sigma12) %*% (ZSig12 %*% 
                                  solve( a = R, b = solve(a = t(R), b = t(ZSig12) %*% (ff - mu)))) )
  
  ## compute the posterior variance of u, the GP at the knot location
  u_var <- Sigma22 + TT + TT %*% solve(a = Sigma22 - TT, b = TT)

  # print("printing min and max ff values at the end of the NR optimization")
  # print(c(min(ff),max(ff)))
  # print("covariance parameters")
  # print(cov_par)
  
  return(list("gp" = ff,"objective_function_values" = obj_fun_vals, "gradient" = grad_psi,
              "u_posterior_mean" = u_mean, "u_posterior_variance" = u_var))
  
}

## test NR algorithm
# cov_par$tau = 0
# ff <- log(c(1,2,3))
# test <- newtrap_sparseGP(start_vals = ff, obj_fun = obj_fun, grad_loglik_fun = grad_loglik_fn, d2log_py_dff = d2log_py_dff,
#                  maxit = 1000, cov_par = cov_par, cov_fun = mv_cov_fun_sqrd_exp, xy = xy, xu = xu, y = y, mu = mu, m = m)

## Test newton-raphson algorithm on simulated sgcp
# source("simulate_sgcp.R")
# set.seed(1308)
# test_sgcp <- simulate_sgcp(lambda_star = 10, gp_par = list("sigma" = 3, "l" = 1, "tau" = sqrt(1e-5)), dbounds = c(0,10), cov_fun = cov_fun_sqrd_exp)
# test_sgcp
# length(test_sgcp)
# 
# plot(x = test_sgcp, y = rep(0, times = length(test_sgcp)))
# hist(test_sgcp, breaks = 100)
# 
# cov_par <- list("sigma" = 3, "l" = 1, "tau" = sqrt(1e-5))
# xy <- seq(from = 0, to = 10 - 0.1, by = 0.2) + 0.1
# xu <- 1:10
# y <- numeric(length(xy))
# for(i in 1:length(xy))
# {
#   y[i] <- sum(test_sgcp > (xy[i] - 0.1) & test_sgcp < (xy[i] + 0.1))
# }
# y
# 
# plot(xy, y)
# 
# start_vals <- rep(log(sum(y)/10), times = length(xy))
# m <- rep(0.2, times = length(xy))
# mu <- rep(log(sum(y)/10), times = length(xy))
# muu <- rep(log(sum(y)/10), times = length(xu))
# 
# 
# system.time(test <- newtrap_sparseGP(start_vals = fstart, obj_fun = obj_fun_pois, 
#                                      grad_loglik_fn = grad_loglik_fn_pois,
#                                      d2log_py_dff = d2log_py_dff_pois,
#                  maxit = 1000, cov_par = cov_par, cov_fun = "sqexp",
#                  xy = as.matrix(xy, ncol = 1), xu = as.matrix(xu, ncol = 1), y = y, mu = mu, muu = muu, m = m))
# 
# test
# test$objective_function_values
# 
# points(x = xy, y = exp(test$gp * 0.2), col = "red", type = "l")

newtrap_sparseGP_update <- function(ff, W, Z, Sigma12, Sigma22, grad_psi, dlog_py_dff, y, mu, ...)
{
  ## ff are the GP values at the observed data locations
  ## W second derivative of p(y|xu, theta) wrt ff
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations 
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21
  ## grad_psi is the gradient of psi = log( p(y|lambda) ) + log( p(ff|theta, x'))
  
  ## Compute Z^(-1) Sigma12
  ZSig12 <- (1/Z) * Sigma12
  
  ## Compute W - Z^(-1)
  WmZ <- W - (1/Z)
  
  ## Compute the Newton-Raphson update 
  
  ## compute -H^(-1) %*% grad_psi
  R <- chol(Sigma22 + t(Sigma12) %*% ZSig12)
  R3 <- chol(Sigma22 + t(Sigma12) %*% (( (Z - 1/W)^(-1) * Sigma12) ) )
  
  # update_complicated <- - 1/WmZ * grad_psi + 
  #   ((1/WmZ) * ZSig12 %*% solve( (Sigma22 + t(Sigma12) %*% ZSig12) + t(ZSig12) %*% ((1/WmZ) * ZSig12) ) ) %*%
  #   (t(ZSig12) %*% ( (1/WmZ) * grad_psi))
  
  # update_complicated <- (- 1/WmZ * grad_psi +
  #   (t(solve( a = t( (Sigma22 + t(Sigma12) %*% ZSig12) + t(ZSig12) %*% ((1/WmZ) * ZSig12) ) ,
  #             b = t((1/WmZ) * ZSig12) )
  #      )
  #    ) %*%
  #   (t(ZSig12) %*% ( (1/WmZ) * grad_psi)))
  
  # AA <- t(ZSig12) %*% ((1/WmZ) * ZSig12)
  # GG <- diag(nrow(Sigma22)) + solve(a = t(R), b = t(solve(a = t(R), b = t(AA))))

  # update_complicated_old2 <- - 1/WmZ * grad_psi +
  #   (t(
  #      solve(a = R, b = solve( a = GG ,
  #             b = solve(a = t(R), b = t((1/WmZ) * ZSig12) ))
  #      ))
  #    ) %*%
  #   (t(ZSig12) %*% ( (1/WmZ) * grad_psi))
  
  ## try new version of update 
  A11 <- (Z / (1 - Z*W)) * dlog_py_dff(ff = ff, y = y, ...)
  A12 <- (1 - Z * W)^(-1) * (ff - mu)
  A13 <- t(solve(a = t(R), b = t((1 / (1 - Z * W)) * Sigma12))) %*%
    solve(a = t(R), b = t(ZSig12) %*% (ff - mu))

  A2 <- t(solve(a = t(R3), b = t(((1 - Z*W)^(-1)) * Sigma12))) %*%
    solve(a = t(R3), b = t(Sigma12) %*% ( (1 / (1 - Z*W)) * grad_psi))

  # A2v2 <- t(solve(a = t(R4), b = t(((1 - Z*W)^(-1)) * Sigma12))) %*%
  #   solve(a = t(R4), b = t(Sigma12) %*% ( (1 / (1 - Z*W)) * grad_psi))

  update_complicated <- A11 - A12 + A13 + A2
  
  # if(any(is.nan(as.numeric(update_complicated))))
  # {
  #   print("gradient of psi is...")
  #   print(grad_psi)
  #   
  #   print("W is...")
  #   print(W)
  #   
  #   print("Z is...")
  #   print(Z)
  #   
  #   print("WmZ is...")
  #   print(WmZ)
  #   
  #   print("AA is...")
  #   print(AA)
  #   
  #   print("GG is...")
  #   print(GG)
  #   # print(update_complicated)
  #   print("ERROR updating ff")
  #   return()
  # }
  
  # nr_update <- ff + update_complicated
  
  # print("old vs new update difference")
  # print(max(abs(update_complicated)))
  # print(max(abs(update_complicated_old2)))
  
  nr_update <- ff + update_complicated
  return(as.numeric(nr_update))
  
  
}

## test the NR update function
# d2log_py_dff <- function(m, ff)
# {
#   return(-m * exp(ff))
# }
# W <- -m * exp(ff)
# newtrap_sparseGP_update(ff = ff, W = W, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, grad_psi = grad_psi)


## Function to compute inverse of the covariance matrix at the observed data locations
# invert_sigma <- function(xy, xu, cov_par, cov_fun, tau = 0)
# {
#   ## xy is a matrix where the rows correspond to observed data locations
#   ## xu is a matrix where the rows correspond to the knot locations
#   ## cov_par is a list of covariance parameters 
#   ## cov_fun is the covariance function
#   ## tau is the potential nugget standard deviations
#   
#   ## create covariances between observed data and knot locations
#   Sigma12 <- matrix(nrow = nrow(xy), ncol = nrow(xu))
#   for(i in 1:nrow(xy))
#   {
#     for(j in 1:nrow(xu))
#     {
#       Sigma12[i,j] <- cov_fun(x1 = xy[i,], x2 = xu[j,], cov_par = cov_par)
#     }
#   }
#   
#   ## Create variance covariancem matrix at knot locations
#   Sigma22 <- make_cov_mat(x = xu, x_pred = numeric(), cov_fun = cov_fun, cov_par = cov_par, tau = tau)
#   
#   U <- Sigma12 %*% solve(a = Sigma22, b = t(Sigma12))
#   
#   ## create the vector of values that get added to the diagonal of Sigma12 %*% Sigma22^-1 %*% Sigma21
#   Z <- cov_par$sigma^2 + tau^2 - diag(U)
#   
#   
#   
# }

## test the gradient function
# y <- c(1,2,3)
# m <- c(1,1,1)
# ff <- rnorm(n = 3)
# mu <- c(0,0,0)
# xu <- matrix(c(0,0.5), ncol = 1)
# xy <- matrix(c(-1,1,2), ncol = 1)
# cov_par <- list("sigma" = 1, "l" = 1)
# 
# Sigma12 <- matrix(nrow = nrow(xy), ncol = nrow(xu))
# for(i in 1:nrow(xy))
# {
#   for(j in 1:nrow(xu))
#   {
#     Sigma12[i,j] <- mv_cov_fun_sqrd_exp(x1 = xy[i,], x2 = xu[j,], cov_par = cov_par)
#   }
# }
# 
# Sigma22 <- make_cov_mat(x = xu, x_pred = numeric(), cov_fun = mv_cov_fun_sqrd_exp, cov_par = cov_par, tau = 0)
# 
# Z <- 1^2 + 0^2 - diag(Sigma12 %*% solve(a = Sigma22, b = t(Sigma12)))
# 
# grad_psi <- grad_loglik_fn(ff = ff, y = y, m = m, mu = mu, Sigma12 = Sigma12, Sigma22 = Sigma22, Z = Z)



