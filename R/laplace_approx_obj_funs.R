## objective functions for the Laplace approximation

## Gaussian objective function
##  note that this does not require the Laplace approximation
#'@export
obj_fun_norm <- function(ff = NA, mu, Z, Sigma12, Sigma22, y, ...)
{
  ## ff are the (maximized) GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21

  ## make sure ff is numeric
  # ff <- as.numeric(ff)
  y <- as.numeric(y)

  ## second derivative of p(y|xu, theta) wrt ff
  args <- list(...)

  ## compute the quadratic form of the GP contribution to the objective function
  ZSig12 <- (1/Z) * Sigma12

  ## compute cholesky decomposition of Sigma22 + Sigma21 %*% Z_inv Sigma12 and the log of the det
  # print(Sigma22)
  # print(t(Sigma12) %*% ZSig12)
  # print(Sigma22 + t(Sigma12) %*% ZSig12)
  # print(head(diag(Sigma22)))

  # print(head(diag(Sigma22 + t(Sigma12) %*% ZSig12)))
  R <- chol(x = Sigma22 + t(Sigma12) %*% ZSig12)
  logdetR <- 2 * sum(log(diag(R)))


  # quad_form_part <- -(1/2) * (t(y - mu) %*% ( (1/Z) * (y - mu))) +
  #   (1/2) * (t(y - mu) %*% ZSig12) %*% solve(a = Sigma22 + t(Sigma12) %*% ZSig12 , b = t((t(y - mu) %*% ZSig12)))
  quad_form_part <- -(1/2) * (t(y - mu) %*% ( (1/Z) * (y - mu))) +
    (1/2) * (t(y - mu) %*% ZSig12) %*% solve(a = R, b = solve(a = t(R), b = t((t(y - mu) %*% ZSig12))))

  ## compute the contribution from the determinents (CORRECT)
  det_part <- -(1/2) * ( sum(log(Z)) - log(det(Sigma22)) +
                             # log(det( solve(Sigma22) + t(Sigma12) %*% ZSig12 ))
                             logdetR
  )

  return(
    quad_form_part + det_part - (length(y) / 2) * log(2*pi)
  )

}

## objective function for a full GP with Gaussian likelihood
#'@export
obj_fun_norm_full <- function(mu, Sigma, y, ...)
{
  return(
    mvtnorm::dmvnorm(x = y, mean = mu, sigma = Sigma, log = TRUE)
  )
}


## Poisson objective function
## objective function for a full GP with poisson likelihood
#'@export
obj_fun_pois_full <- function(ff, mu, Sigma, y, ...)
{
  ## ff are the (maximized) GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma is the variance covariance matrix at the observed data locations

  ## make sure ff is numeric
  ff <- as.numeric(ff)

  ## second derivative of p(y|xu, theta) wrt ff
  args <- list(...)
  m <- args$m
  W <- -m * exp(ff)
  grad_log_py_ff <- dlog_py_dff_pois(ff = ff, y = y, ...)

  ## log(p(y|xu, theta))
  log_py <- sum(y * log(m) - lfactorial(y) - m * exp(ff) + y * ff)

  ## cholesky decomposition of (I + sqrt(-W) %*% Sigma %*% sqrt(-W))
  L <- t(chol(diag(length(y)) + diag(sqrt(-W)) %*% Sigma %*% diag(sqrt(-W))))

  ## calculate (-W) %*% f + (d/df) log(p(y|f))
  b <- (-diag(W)) %*% (ff - mu) + grad_log_py_ff

  ## calculate b - ( sqrt(-W) %*% t(L) )_inverse %*% (L_inverse sqrt(-W) %*% Sigma b)
  a <- b - sqrt(-diag(W)) %*% solve(a = t(L) ,
                 b = solve(a = L, b = sqrt(-diag(W)) %*% Sigma %*% b))

  return(
    ## pulled straight from page 46 of R & W
    -(1/2) * (t(a) %*% (ff - mu)) + log_py - sum(log(diag(L) ) )
  )

}

## create the the objective function log q(y|theta, xu, ff_hat)
## This function is perfect
## objective function for a full GP with poisson likelihood
#'@export
obj_fun_pois <- function(ff, mu, Z, Sigma12, Sigma22, y, ...)
{
  ## ff are the (maximized) GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21

  ## make sure ff is numeric
  ff <- as.numeric(ff)

  ## second derivative of p(y|xu, theta) wrt ff
  args <- list(...)
  m <- args$m
  W <- -m * exp(ff)
  Z2 <- 1 + sqrt(-W) * Z * sqrt(-W)


  ## log(p(y|xu, theta))
  log_py <- sum(y * log(m) - lfactorial(y) - m * exp(ff) + y * ff)

  ## compute the quadratic form of the GP contribution to the objective function
  ZSig12 <- (1/Z) * Sigma12

  ## compute cholesky decomposition of Sigma22 + Sigma21 %*% Z_inv Sigma12 and the log of the det
  # print(Sigma22)
  # print(t(Sigma12) %*% ZSig12)
  # print(Sigma22 + t(Sigma12) %*% ZSig12)
  # print(head(diag(Sigma22)))

  # print(head(diag(Sigma22 + t(Sigma12) %*% ZSig12)))
  R <- chol(x = Sigma22 + t(Sigma12) %*% ZSig12)
  logdetR <- 2 * sum(log(diag(R)))

  ## NEW
  R2 <- chol(x = Sigma22 + t(sqrt(-W) * Sigma12) %*% ((1 / Z2) * (sqrt(-W) * Sigma12)))
  logdetR2 <- 2 * sum(log(diag(R2)))

  ## take cholesky decomp of sigma22
  R_Sigma22 <- chol(x = Sigma22)

  ## compute cholesky decomposition of Sigma22 + t(Sigma12) %*% ZSig12 - t(ZSig12) %*% ((1/(-W + 1/Z) * ZSig12) ) and the log of the det
  ## See Vanhatalo et al. 2010

  # quad_form_part <- -(1/2) * (t(ff - mu) %*% ( (1/Z) * (ff - mu))) +
  #   (1/2) * (t(ff - mu) %*% ZSig12) %*% solve(a = Sigma22 + t(Sigma12) %*% ZSig12 , b = t((t(ff - mu) %*% ZSig12)))
  # quad_form_part <- -(1/2) * (t(ff - mu) %*% ( (1/Z) * (ff - mu))) +
  #   (1/2) * (t(ff - mu) %*% ZSig12) %*% solve(a = R, b = solve(a = t(R), b = t((t(ff - mu) %*% ZSig12))))
  quad_form_part <- -(1/2) * (t(ff - mu) %*% ( (1/Z) * (ff - mu))) +
    (1/2) * (t(solve(a = t(R), b = t((t(ff - mu) %*% ZSig12))))) %*% solve(a = t(R), b = t((t(ff - mu) %*% ZSig12)))


  ## compute the contribution from the determinents
  ## NEW
  det_part_1 <- -(1/2) * ( - 2*sum(log(diag(R_Sigma22))) +
                             # log(det( solve(Sigma22) + t(Sigma12) %*% ZSig12 ))
                             logdetR2
  )
  det_part_2 <- -(1/2) * sum(log(Z2))


  return(
    quad_form_part + log_py + det_part_1 + det_part_2
  )
}

## logistic function
#'@export
my_logistic <- function(x)
{
    return(
      1 / (1 + exp(-x))
    )
}

## Bernoulli objective function
## create the the objective function log q(y|theta, xu, ff_hat)
## This function is likely okay
#'@export
obj_fun_bern <- function(ff, mu, Z, Sigma12, Sigma22, y, ...)
{
  ## ff are the (maximized) GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21

  ## make sure ff is numeric
  ff <- as.numeric(ff)

  # sigma <- cov_par$sigma

  ## second derivative of p(y|xu, theta) wrt ff

  ## compute probability vector pi(ff)
  # pi_ff <- (1 - 1e-8) * my_logistic(x = ff) + (1e-8) * 0.5
  pi_ff <- my_logistic(x = ff)
  log_pi_ff <- plogis(q = ff, log.p = TRUE)
  log_1minus_pi_ff <- plogis(q = -ff, log.p = TRUE)


  ## compute first derivative vector of link wrt ff
  dpi_dff <- pi_ff * (1 - pi_ff)
  # log_dpi_dff <- log_pi_ff + log_1minus_pi_ff

  ## compute second derivative vector of link wrt ff
  d2pi_dff2 <- dpi_dff * (1 - 2 * pi_ff)
  # d2pi_dff2 <- dpi_dff - 2 * exp(log_pi_ff + log_dpi_dff)

  ## get matrix W
  args <- list(...)
  # W2 <- (y / pi_ff - 1 / (1 - pi_ff) + y / (1 - pi_ff)) * d2pi_dff2 +
  #   (dpi_dff)^2 * (-y / (pi_ff^2) - 1 / ((1 - pi_ff)^2) - y / ((1 - pi_ff)^2) )
  W <- (1 - 2*pi_ff) * (y - pi_ff) -
    (y*(1 - pi_ff)^2 + pi_ff^2 + y * pi_ff^2)

  if(any(is.nan(W)) || (max(abs(ff)) > 1e100 ))
  {
    print(c(min(ff), max(ff)))
    # W[is.nan(W)] <- 0
    print("second derivatives of the likelihood wrt f encountered an error")
  }

  if(any(W > 0))
  {
    print("There is some W > 0")
  }

  # print(summary(pi_ff))
  # print(summary(dpi_ff))
  # print(summary(pi_ff))

  # print(c(min(pi_ff), max(pi_ff)))

  Z2 <- (1 + sqrt(-W) * Z * sqrt(-W))

  ## log(p(y|xu, theta))
  log_py <- sum(
    y * log_pi_ff + (1 - y) * log_1minus_pi_ff
  )

  ## compute the quadratic form of the GP contribution to the objective function
  ZSig12 <- ((1/Z) * Sigma12 )

  ## compute cholesky decomposition of Sigma22 + Sigma21 %*% Z_inv Sigma12 and the log of the det
  R <- chol(x = Sigma22 + t(Sigma12) %*% ZSig12)
  logdetR <- 2 * sum(log(diag(R)))

  ## NEW
  R2 <- chol(x = Sigma22 + t(sqrt(-W) * Sigma12) %*% ((1 / Z2) * (sqrt(-W) * Sigma12)))
  logdetR2 <- 2 * sum(log(diag(R2)))

  ## cholesky decomp of Sigma22
  R_sigma22 <- chol(x = Sigma22)


  ## compute cholesky decomposition of Sigma22 + t(Sigma12) %*% ZSig12 - t(ZSig12) %*% ((1/(-W + 1/Z) * ZSig12) ) and the log of the det

  ## DEPRECATED
  # V <- diag(nrow(Sigma22)) - solve(a = t(R), b = t(ZSig12)) %*% (( (1/(-W + 1/Z)) * t(solve(a = t(R), b = t(ZSig12))) ) )
  # cholV <- chol(x = V)
  # old_logdetR2 <- logdetR + 2 * sum(log(diag(cholV)))
  # old_R2 <- chol(x = Sigma22 + t(Sigma12) %*% ZSig12 - t(ZSig12) %*% (( (1/(-W + 1/Z)) * ZSig12) ) )
  # old_logdetR2 <- 2 * sum(log(diag(old_R2)))

  # quad_form_part <- -(1/2) * (t(ff - mu) %*% ( (1/Z) * (ff - mu))) +
  #   (1/2) * (t(ff - mu) %*% ZSig12) %*% solve(a = Sigma22 + t(Sigma12) %*% ZSig12 , b = t((t(ff - mu) %*% ZSig12)))
  # quad_form_part <- -(1/2) * (t(ff - mu) %*% ( (1/Z) * (ff - mu))) +
  #   (1/2) * (t(ff - mu) %*% ZSig12) %*% solve(a = R, b = solve(a = t(R), b = t((t(ff - mu) %*% ZSig12))))
  quad_form_part <- -(1/2) * (t(ff - mu) %*% ( (1/Z) * (ff - mu))) +
    (1/2) * (t(solve(a = t(R), b = t((t(ff - mu) %*% ZSig12))))) %*% solve(a = t(R), b = t((t(ff - mu) %*% ZSig12)))

  ## compute the contribution from the determinents
  ## This part can be seen in Vanhatalo et al. 2010
  # old_det_part_1 <- -(1/2) * (sum(log(Z)) - log(det(Sigma22)) +
  #                            logdetR
  # )

  ## NEW
  # det_part_1 <- -(1/2) * ( - log(det(Sigma22)) +
  #                            # log(det( solve(Sigma22) + t(Sigma12) %*% ZSig12 ))
  #                            logdetR2
  # )
  det_part_1 <- -(1/2) * ( - 2*sum(log(diag(R_sigma22))) +
                             # log(det( solve(Sigma22) + t(Sigma12) %*% ZSig12 ))
                             logdetR2
  )
  det_part_2 <- -(1/2) * sum(log(Z2))




  # old_det_part_2 <- -(1/2) * ( sum(log(1/Z - W)) -
  #                            # log(det(Sigma22 + t(Sigma12) %*% ZSig12)) +
  #                            logdetR +
  #                            old_logdetR2)

  # if(any(W == 0))
  # {
  #   print("objective function value is...")
  #   print(quad_form_part + log_py + det_part_1 + det_part_2)
  #
  #   print("min and max probabilities are...")
  #   print(c(min(pi_ff), max(pi_ff)))
  #
  #   print("min and max GP values are...")
  #   print(c(min(ff),max(ff)))
  # }
  # if(is.nan(quad_form_part) || quad_form_part == Inf || quad_form_part == -Inf)
  # {
  #   print("quad_form_part")
  # }
  # if(is.nan(log_py) || log_py == Inf || log_py == -Inf)
  # {
  #   print("log_py")
  #   print(log_py)
  # }
  # if(is.nan(det_part_1) || det_part_1 == Inf || det_part_1 == -Inf)
  # {
  #   print("det_part1")
  # }
  # if(is.nan(det_part_2) || det_part_2 == Inf || det_part_2 == -Inf)
  # {
  #   print("det_part2")
  # }
  return(
    quad_form_part + log_py + det_part_1 + det_part_2
  )

}

## create the the objective function for FULL GP with Bernoulli data
## This function is likely okay
#'@export
obj_fun_bern_full <- function(ff, mu, Sigma, y, ...)
{
  ## ff are the (maximized) GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma variance covariance matrix at observed data locations

  ## make sure ff is numeric
  ff <- as.numeric(ff)

  ## second derivative of p(y|xu, theta) wrt ff

  ## compute probability vector pi(ff)
  pi_ff <- my_logistic(x = ff)
  log_pi_ff <- plogis(q = ff, log.p = TRUE)
  log_1minus_pi_ff <- plogis(q = -ff, log.p = TRUE)

  ## compute first derivative vector of link wrt ff
  dpi_dff <- pi_ff * (1 - pi_ff)

  ## compute second derivative vector of link wrt ff
  d2pi_dff2 <- dpi_dff * (1 - 2 * pi_ff)

  grad_log_py_ff <- dlog_py_dff_bern(ff = ff, y = y, ...)

  ## get matrix W
  args <- list(...)
  W <- (1 - 2*pi_ff) * (y - pi_ff) -
    (y*(1 - pi_ff)^2 + pi_ff^2 + y * pi_ff^2)

  ## log(p(y|xu, theta))
  log_py <- sum(
    y * log_pi_ff + (1 - y) * log_1minus_pi_ff
  )

  ## cholesky decomposition of (I + sqrt(-W) %*% Sigma %*% sqrt(-W))
  L <- t(chol(diag(length(y)) + diag(sqrt(-W)) %*% Sigma %*% diag(sqrt(-W))))

  ## calculate (-W) %*% f + (d/df) log(p(y|f))
  b <- (-diag(W)) %*% (ff - mu) + grad_log_py_ff

  ## calculate b - ( sqrt(-W) %*% t(L) )_inverse %*% (L_inverse sqrt(-W) %*% Sigma b)
  a <- b - sqrt(-diag(W)) %*% solve(a = t(L) ,
                                    b = solve(a = L, b = sqrt(-diag(W)) %*% Sigma %*% b))

  return(
    ## pulled straight from page 46 of R & W
    -(1/2) * (t(a) %*% (ff - mu)) + log_py - sum(log(diag(L) ) )
  )

}

## objective function wrapper to be passed to optim
#'@export
obj_fun_wrapper <- function(cov_par, ...)
{
  ## cov_par is a vector of the parameters being optimized over
  ## other_cov_par is a vector of covariance parameters not being optimized over
  ## ffstart are the starting values for the GP at the observed data locations
  ## y is the vector of observed data
  ## mu is a vector of GP means of the same length as ff
  ## muu is a vector of GP means at knot locations of length nrow(xu)
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21
  ## maxit
  ## tol_nr
  ## grad_loglik_fn
  ## d2log_py_dff
  ## obj_fun
  ## transform
  ## dcov_fun_dtheta
  ## num_cov_par is the number of elements in cov_par that correspond to covariance parameters that are not knots
  ##  these parameters must always come first in cov_par
  args <- list(...)

  # ff <- args$ff
  # y <- args$y
  # mu <- args$mu
  # muu <- args$muu
  # Sigma12 <- args$Sigma12

  ## store back transformed transformed parameter values
  trans_cov_par <- cov_par
  if(args$transform == TRUE)
  {
    for(p in 1:(args$num_cov_par))
    {
      ## this is sort of a hacky line that just transforms parameters back to the right scales
      trans_cov_par[p] <- dcov_fun_dtheta[[p]](x1 = 0, x2 = 2,
                                           cov_par = as.list(abs(cov_par[1:args$num_cov_par])),
                                           transform = args$transform)$trans_fun(cov_par[p])
    }
  }

  print(trans_cov_par)

  ## find the value of the GP that maximizes log(p(y|f)) + log(p(f|theta, xu))
  if(is.function(args$dcov_fun_dknot))
  {
    xu <- matrix(nrow = nrow(args$xustart), ncol = ncol(args$xustart), data = cov_par[(args$num_cov_par + 1):length(cov_par)], byrow = TRUE)
  }
  if(!is.function(args$dcov_fun_dknot))
  {
    xu <- args$xustart
  }
  nr_step <- newtrap_sparseGP(start_vals = args$ffstart,
                              cov_par = as.list(c(trans_cov_par[1:args$num_cov_par], args$other_cov_par)),
                              xu = xu,
                              ...)

  return(nr_step$objective_function_values[length(nr_step$objective_function_values)])

}

# ## Test objective function to be passed to optim
# set.seed(1308)
# cov_par <- list("sigma" = 3, "l" = 1, "tau" = sqrt(1e-5))
# xy <- matrix(
#   c(rep(seq(from = 0.01, to = 9.98, by = area), times = 1)),
#   ncol = 1, byrow = FALSE)
# # xu <- matrix(c(rep(seq(from = 0, to = 10, by = 1), times = 1)),
# #              ncol = 1, byrow = FALSE)
# mu <- rep(0, times = nrow(xy))
# # m <- rep(area, times = nrow(xy))
# cov_fun <- mv_cov_fun_sqrd_exp
#
# # Sigma22 <- make_cov_mat(x = xu, x_pred = numeric(), cov_fun = cov_fun, cov_par = cov_par)
# Sigma <- make_cov_mat(x = xy, x_pred = numeric(), cov_fun = cov_fun, cov_par = cov_par, tau = cov_par$tau)
# # log_lambda_u <- as.numeric(mvtnorm::rmvnorm(n = 1, sigma = Sigma22))
# # cond_mean_y <- as.numeric(normal_cond_mean(y = log_lambda_u, x = xu, x_pred = xy,
# #                                 mu = rep(0, length = nrow(rbind(xu, xy))), sigma = Sigma))
#
# # cond_var_y <- diag(normal_cond_var( x = xu, x_pred = xy, sigma = Sigma))
#
# log_lambda_y <- mvtnorm::rmvnorm(n = 1, sigma = Sigma)
#
# plot(xy, log_lambda_y, type = "l")
# y <- rpois(n = length(log_lambda_y), lambda = exp(log_lambda_y) * rep(area, times = nrow(xy)))
#
# dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl)
# dcov_fun_dknot <- dsqexp_dx2
#
# xu <- matrix(nrow = length(c(seq(from = 1, to = 8, by = 1))), ncol = 1, data = c(seq(from = 1, to = 8, by = 1)))
# muu <- rep(0, times = nrow(xu))
#
# cov_par_start <- list("sigma" = log(4), "l" = log(1.28906354), "tau" = sqrt(1e-5))
# fstart <- rep(log(max(y)) - log(area), times = nrow(xy))
#
# cov_par_vec <- unlist(cov_par_start)[1:2]
# # cov_par_vec <- c(cov_par_vec, t(xu))
# other_cov_par <- c("tau" = sqrt(1e-5))
# obj_fun_wrapper(cov_par = cov_par_vec,
#                 other_cov_par = other_cov_par,
#                 num_cov_par = 2,
#                 cov_fun = mv_cov_fun_sqrd_exp,
#                 dcov_fun_dtheta = dcov_fun_dtheta,
#                 dcov_fun_dknot = NA,
#                 xustart = xu,
#                 xy = xy,
#                 y = y,
#                 ffstart = fstart,
#                 grad_loglik_fn = grad_loglik_fn_pois,
#                 dlog_py_dff = dlog_py_dff_pois,
#                 d2log_py_dff = d2log_py_dff_pois,
#                 mu = mu,
#                 muu = muu,
#                 transform = TRUE,
#                 obj_fun = obj_fun_pois,
#                 maxit = 1000,
#                 tol = 1e-5,
#                 verbose = TRUE,
#                 m = rep(area, times = nrow(xy)))
