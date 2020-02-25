## My Newton-Raphson algorithm for finding the mode
##  of the Laplace approximation to a FULL GP model with non-Gaussian likelihood
#source("gp_functions.R")
# Rcpp::sourceCpp(file = "covariance_functions.cpp")
#'@export
newtrapGP <- function(start_vals,
                             obj_fun,
                             grad_loglik_fn,
                             dlog_py_dff,
                             d2log_py_dff,
                             maxit = 10000,
                             tol = 1e-6,
                             cov_par,
                             cov_fun,
                             xy,
                             y,
                             mu,
                             delta = 1e-6,
                             ...)
{
  ## start_vals is a vector of starting values of the GP corresponding to the observed data locations
  ## maxit is the maximum number of iterations of the algorithm
  ## tol is the absolute tolerance level which dictates
  ##      how small a difference the algorithm is allowed to make before stopping
  ##      to the objective function as well as how close the gradient is allowed to be to zero
  ## obj_fun is the objective function
  ## grad_loglik_fun is a function that computes the gradient of
  ##      psi = log(p(y|lambda)) + log(p(ff|theta)) wrt ff
  ## d2log_py_dff is a function that computes the vector of second derivatives
  ##      of log(p(y|lambda)) wrt ff
  ## dlog_py_dff is a function that computes the vector of second derivatives
  ##      of log(p(y|lambda)) wrt ff
  ## cov_par is a list of the covariance function parameters
  ## cov_fun is a string specifying either "sqexp" or "exp"
  ## xy are the observed data locations
  ## xu are the unobserved knot locations
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## muu is the mean of the GP at each of the knot locations
  ## ... are argument to be passed to the objective function and the gradient of the log-likelihood

  ## create Sigma
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
    Sigma <- make_cov_mat_ardC(x = xy, x_pred = matrix(),
                           cov_fun = cov_fun,
                           cov_par = cov_par,
                           delta = delta,
                           lnames = lnames)

  }
  else{
    Sigma <- make_cov_matC(x = xy, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)

  }
  ## initialize vector of objective function values
  ff <- start_vals
  obj_fun_vals <-  numeric()
  obj_fun_vals[1] <- obj_fun(ff = ff, Sigma = Sigma, y = y, mu = mu, ...)

  ## run NR updates until convergence criteria is met or the maxit is reached
  iter <- 1

  ## update the iteration counter
  iter <- iter + 1

  ## update ff
  W <- as.numeric(d2log_py_dff(ff = ff, y = y, ...))
  grad_log_py_ff <- dlog_py_dff(ff = ff, y = y, ...)
  grad_psi <- grad_loglik_fn(ff = ff, y = y, mu = mu, Sigma = Sigma, ...)
  ff <- newtrapGP_update(ff = ff, mu = mu, W = W, Z = Z, Sigma = Sigma, grad_log_py_ff = grad_log_py_ff, ...)

  ## compute objective function
  obj_fun_vals[iter] <- obj_fun(ff = ff, Sigma = Sigma, y = y, mu = mu, ...)

  while(iter < maxit && (abs(obj_fun_vals[iter] - obj_fun_vals[iter - 1]) > tol || any(abs(grad_psi) > tol)))
  {
    ## update the iteration counter
    iter <- iter + 1

    ## update ff
    grad_psi <- grad_loglik_fn(ff = ff, y = y, mu = mu, Sigma = Sigma, ...)
    grad_log_py_ff <- dlog_py_dff(ff = ff, y = y, ...)
    W <- d2log_py_dff(ff = ff, y = y,  ...)
    ff <- newtrapGP_update(ff = ff, mu = mu,  W = W, Sigma = Sigma, grad_log_py_ff = grad_log_py_ff, ...)

    ## compute objective function
    obj_fun_vals[iter] <- obj_fun(ff = ff, Sigma = Sigma, y = y, mu = mu, ...)
  }

  ## return matrices used for prediction from this Laplace approximation
  L <- t(chol(diag(length(ff)) + sqrt(-diag(W)) %*% Sigma %*% sqrt(-diag(W)) ))


  return(list("gp" = ff,"objective_function_values" = obj_fun_vals,
              "gradient" = grad_psi,
              "grad_log_py_ff" = grad_log_py_ff,
              "L" = L, "W" = W))

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

## computations straight out of pg 46 of R & W
#'@export
newtrapGP_update <- function(ff, W, Sigma, grad_log_py_ff, y, mu, ...)
{
  ## ff are the GP values at the observed data locations
  ## W second derivative of p(y|xu, theta) wrt ff
  ## Sigma is the variance covariance matrix at the observed data locations
  ## grad_log_py_ff is the gradient of the likelihood wrt f
  ## y is the observed response vector
  ## ... are arguments to be passed to the dlog_py_dff

  ## Compute the Newton-Raphson update
  ff <- as.numeric(ff)

  ## compute cholesky decomposition of I + sqrt(-diag(W)) %*% Sigma %*% sqrt(-diag(W))
  L <- t(chol(diag(length(ff)) + sqrt(-diag(W)) %*% Sigma %*% sqrt(-diag(W))))

  ## calculate (-W) %*% f + (d/df) log(p(y|f))
  b <- (-diag(W)) %*% (ff - mu) + grad_log_py_ff

  ## calculate b - ( sqrt(-W) %*% t(L) )_inverse %*% (L_inverse sqrt(-W) %*% Sigma b)
  a <- b - sqrt(-diag(W)) %*% solve(a = t(L) ,
                                    b = solve(a = L, b = sqrt(-diag(W)) %*% Sigma %*% b))

  nr_update <- Sigma %*% a + mu

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

