## knot proposal functions
#'@export
knot_prop_var <- function(laplace_opt, ...)
{
  ## laplace_opt is a list of output returned from laplace_grad_ascent
  ## ... should contain the function to get the posterior predictive variances at the data locations
  ##        that function is called predict_laplace()

  pred_vars <- predict_laplace(u_mean = laplace_opt$u_mean,
                              u_var = laplace_opt$u_var,
                              xu = laplace_opt$xu,
                              x_pred = laplace_opt$xy,
                              cov_fun = laplace_opt$cov_fun,
                              cov_par = laplace_opt$cov_par,
                              mu = laplace_opt$mu,
                              muu = laplace_opt$muu)$pred_var
  knot_prop <- matrix(nrow = 1, laplace_opt$xy[which.max(as.numeric(diag(pred_vars))),])
  return(knot_prop)
}

## knot prop pred error
#'@export
knot_prop_error <- function(laplace_opt, ...)
{
  ## laplace_opt is a list of output returned from laplace_grad_ascent
  ## ... should contain the function to get the posterior predictive variances at the data locations
  ##        that function is called predict_laplace()

  preds <- predict_laplace(u_mean = laplace_opt$u_mean,
                               u_var = laplace_opt$u_var,
                               xu = laplace_opt$xu,
                               x_pred = laplace_opt$xy,
                               cov_fun = laplace_opt$cov_fun,
                               cov_par = laplace_opt$cov_par,
                               mu = laplace_opt$mu,
                               muu = laplace_opt$muu)$pred_mean
  error <- abs(exp(laplace_opt$fmax) - exp(preds))

  knot_prop <- matrix(nrow = 1, laplace_opt$xy[which.max(error),])
  return(knot_prop)
}

## knot prop EGO
#'@export
knot_prop_ego <- function(laplace_opt,
                          # obj_fun,
                          # grad_loglik_fn,
                          # d2log_py_dff,
                          # maxit_nr = 2000,
                          # tol_nr = 1e-5,
                          # y,
                          opt = NA, ...)
{
  ## laplace_opt is a list of output returned from laplace_grad_ascent
  ## obj_fun is the objective function (log(p(y|f) + log(p(f|theta, xu)))
  ## d2log_py_dff is a function that computes the vector of second derivatives
  ##      of log(p(y|lambda)) wrt ff
  ## grad_loglik_fun is a function that computes the gradient of
  ##      psi = log(p(y|lambda)) + log(p(ff|theta)) wrt ff
  ## ... should contain the covariance function to be used to make proposals
  ##      call this function ego_cov_fun(x1, x2, cov_par)
  ## ... Must also contain the argument TTmin denoting the minimum number of allowed objective function evaluations
  ## ... Must also contain the argument TTmin denoting the maximum number of allowed objective function evaluations
  ## ... Must also contain the argument ego_cov_par: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_cov_fun: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_dcov_fun_dtheta: a list of functions which corresponds
  ##        (in order) to the covariance parameters in
  ##        cov_par. These functions compute the derivatives of the
  ##        covariance function wrt the particular covariance parameter.
  ## ... Must also contain the argument predict_laplace: to get predictive variances for the first prediction

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 0.8,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = 1e-1, "maxit_nr" = 1000, "delta" = 1e-6,
                    "tol_knot" = 1e-1, "maxknot" = 30, "tol_nr" = 1e-4, "TTmin" = 10, "TTmax" = 40,
                    "ego_cov_par" = list("sigma" = 3, "l" = 1, "tau" = 0), "ego_cov_fun" = "exp",
                    "ego_cov_fun_dtheta" = list("sigma" = dexp_dsigma,
                                                "l" = dexp_dl))

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: invalid or unnecessary argument to opt(). Proceeding with fingers crossed.")
        next
      }

      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  optim_method = opt_master$optim_method
  optim_par = list("decay" = opt_master$decay,
                   "epsilon" = opt_master$epsilon,
                   "eta" = opt_master$eta,
                   "learn_rate" = opt_master$learn_rate)
  maxit = opt_master$maxit
  maxit_nr <- opt_master$maxit_nr
  obj_tol = opt_master$obj_tol
  grad_tol = opt_master$grad_tol
  delta <- opt_master$delta
  tol_knot <- opt_master$tol_knot
  maxknot <- opt_master$maxknot
  tol_nr <- opt_master$tol_nr
  TTmin <- opt_master$TTmin
  TTmax <- opt_master$TTmax
  ego_cov_par <- opt_master$ego_cov_par
  ego_cov_fun <- opt_master$ego_cov_fun
  ego_dcov_fun_dtheta <- opt_master$ego_cov_fun_dtheta


  args <- list(...)
  obj_fun <- args$obj_fun
  d2log_py_dff <- args$d2log_py_dff
  dlog_py_dff <- args$dlog_py_dff
  grad_loglik_fun <- args$grad_loglik_fun
  predict_laplace <- args$predict_laplace
  cov_fun <- args$cov_fun

  ## store current knot locations
  xu <- laplace_opt$xu
  cov_par <- laplace_opt$cov_par
  xy <- laplace_opt$xy

  ## store objective function values
  # obj_fun_vals <- rep(laplace_opt$obj_fun[length(laplace_opt$obj_fun)], times = nrow(xu))
  obj_fun_vals <- numeric()

  ## first set of proposals for the new knot
  # pred_vars <- as.numeric(diag(predict_laplace(u_mean = laplace_opt$u_mean,
  #                          u_var = laplace_opt$u_var,
  #                          xu = laplace_opt$xu,
  #                          x_pred = laplace_opt$xy,
  #                          cov_fun = laplace_opt$cov_fun,
  #                          cov_par = laplace_opt$cov_par,
  #                          mu = laplace_opt$mu,
  #                          muu = laplace_opt$muu)$pred_var))


  pseudo_prop <- matrix(laplace_opt$xy[sample.int(n = nrow(laplace_opt$xy), size = TTmin, replace = FALSE),],
                        ncol = ncol(laplace_opt$xy), nrow = TTmin)

  ## meta model x values
  for(i in 1:nrow(pseudo_prop))
  {
    pseudo_xu <- rbind(xu, pseudo_prop[i,])
    # print(pseudo_xu)

    ## meta model y values
    nr_opt <- try(newtrap_sparseGP(start_vals = laplace_opt$fmax,
                               # obj_fun = obj_fun,
                               # grad_loglik_fn = grad_loglik_fn,
                               # dlog_py_dff = dlog_py_dff,
                               # d2log_py_dff = d2log_py_dff,
                               # maxit = maxit_nr,
                               # tol = tol_nr,
                               cov_par = laplace_opt$cov_par,
                               # cov_fun = cov_fun,
                               xy = xy,
                               xu = pseudo_xu,
                               # y = y,
                               mu = laplace_opt$mu,
                               muu = c(laplace_opt$muu, laplace_opt$muu[1]), delta = delta, ...))
    if(class(nr_opt) == "try-error")
    {
      while(class(nr_opt) == "try-error")
      {
        pseudo_prop[i,] <- laplace_opt$xy[sample.int(n = nrow(laplace_opt$xy), size = 1, replace = FALSE),] +
          rnorm(n = length(pseudo_prop[i,]), mean = 0, sd = laplace_opt$cov_par$l / 100)
        pseudo_xu <- rbind(xu, pseudo_prop[i,])
        # print(pseudo_xu)

        ## meta model y values
        nr_opt <- try(newtrap_sparseGP(start_vals = laplace_opt$fmax,
                                       # obj_fun = obj_fun,
                                       # grad_loglik_fn = grad_loglik_fn,
                                       # dlog_py_dff = dlog_py_dff,
                                       # d2log_py_dff = d2log_py_dff,
                                       # maxit = maxit_nr,
                                       # tol = tol_nr,
                                       cov_par = laplace_opt$cov_par,
                                       # cov_fun = cov_fun,
                                       xy = xy,
                                       xu = pseudo_xu,
                                       # y = y,
                                       mu = laplace_opt$mu,
                                       muu = c(laplace_opt$muu, laplace_opt$muu[1]), delta = delta, ...))
      }
    }

    obj_fun_vals <- c(obj_fun_vals, nr_opt$objective_function_values[length(nr_opt$objective_function_values)])
  }
  # obj_fun_x <- rbind(xu, pseudo_prop)
  obj_fun_x <- rbind(pseudo_prop)


  ## get meta model parameters

  ## meta model GP negative log likelihood function
  meta_mod_obj_fun <- function(cov_par, ...)
  {
    # y, x, cov_fun, mu
    args <- list(...)
    x <- args$x
    obj_fun_vals <- args$obj_fun_vals
    cov_fun <- args$cov_fun
    mu <- args$mu
    cov_par_names <- args$cov_par_names
    tau <- args$tau
    delta <- args$delta

    cov_par <- as.list(exp(cov_par))
    names(cov_par) <- cov_par_names
    cov_par$tau <- tau

    # print(cov_par)
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x), sep = "")
      Sig <- make_cov_mat_ardC(x = x, x_pred = matrix(),
                               cov_fun = cov_fun, cov_par = cov_par,
                               delta = delta, lnames = lnames)

    }
    else{
      Sig <- make_cov_matC(x = x, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)

    }
    # print(Sig)

    L <- t(chol(Sig))
    mu <- rep(mu,times = length(obj_fun_vals))

    temp <- -1 * (-1 * sum(log(diag(L))) - 1/2 * t(solve(a = L, b = obj_fun_vals - mu)) %*% solve(a = L, b = obj_fun_vals - mu))

    return(temp)
  }

  ## function to create d(Sigma)/dtheta
  dsig_dtheta <- function(cov_par, par_name, dcov_fun_dtheta, x, transform = TRUE)
  {
    mat <- matrix(nrow = nrow(x), ncol = nrow(x))
    for(i in 1:nrow(mat))
    {
      for(j in 1:ncol(mat))
      {
        temp <- eval(
          substitute(expr = dcov_fun_dtheta$par(x1 = a, x2 = b, cov_par = c, transform = d),
                     env = list("a" = x[i,], "b" = x[j,], "c" = cov_par, "par" = par_name, "d" = transform))
        )
        mat[i,j] <- temp$derivative
      }
    }
    return(mat)
  }

  ## meta model gradient function
  meta_mod_grad <- function(cov_par, ...)
  {
    # y, x, cov_fun, mu, dcov_fun_dtheta, transform = TRUE
    args <- list(...)
    x <- args$x
    obj_fun_vals <- args$obj_fun_vals
    cov_fun <- args$cov_fun
    dcov_fun_dtheta <- args$dcov_fun_dtheta
    transform <- args$transform
    cov_par <- as.list(exp(cov_par))
    cov_par_names <- args$cov_par_names
    names(cov_par) <- cov_par_names
    mu <- args$mu
    mu <- rep(mu, times = nrow(x))
    tau <- args$tau
    cov_par$tau <- tau
    delta <- args$delta




    ## store gradient values
    grad <- numeric(length(dcov_fun_dtheta))

    ## create current covariancne matrix
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x), sep = "")
      Sig <- make_cov_mat_ardC(x = x, x_pred = matrix(), cov_fun = cov_fun,
                               cov_par = cov_par,
                               delta = delta,
                               lnames = lnames)

    }
    else{
      Sig <- make_cov_matC(x = x, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)

    }
    L <- t(chol(Sig))

    ## loop through each covariance parameter
    for(p in 1:length(dcov_fun_dtheta))
    {
      par_name <- names(cov_par)[p]

      ## next compute dSigma/dtheta_j
      dSigma_dtheta <- dsig_dtheta(par_name = par_name,
                                       cov_par = cov_par,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       x = x,
                                       transform = transform)
      grad[p] <- -(
        -1/2 * sum(diag( solve(a = t(L), b = solve(a = L, b = dSigma_dtheta)) )) +
          1/2 * t(solve(a = L, b = obj_fun_vals - mu)) %*% solve(a = L, b = dSigma_dtheta) %*% solve(a = t(L), b = solve(a = L, b = obj_fun_vals - mu))
      )
    }

    return(grad)
  }

  ## optimize meta model covariance parameters
  temp_opt <- try(optim(par = log(as.numeric(ego_cov_par))[1:2], fn = meta_mod_obj_fun, gr = meta_mod_grad, method = "BFGS",
                    obj_fun_vals = obj_fun_vals, x = obj_fun_x, cov_fun = ego_cov_fun,
                    mu = obj_fun_vals[1], dcov_fun_dtheta = ego_dcov_fun_dtheta,
                    transform = TRUE, cov_par_names = names(ego_cov_par)[1:2], tau = ego_cov_par$tau, delta = delta))

  ## If optim hits values of sigma that cause infinite variance or zero length scale, then set parameter values
  ##  to something semi reasonable
  if(class(temp_opt) == "try-error")
  {
    print("Error: Optimizing meta GP was unsuccessful - likely due to extreme parameter values. Setting covariance parameters manually.")
    temp_opt <- list()
    temp_opt$par <- c(log(max(obj_fun_vals) - min(obj_fun_vals)), log( (max(obj_fun_x) - min(obj_fun_x))/1000 ))
  }

  cov_par_names <- names(ego_cov_par)
  ego_cov_par <- as.list(c(exp(temp_opt$par), ego_cov_par$tau))
  names(ego_cov_par) <- cov_par_names

  K12 <- make_cov_matC(x = laplace_opt$xy, x_pred = obj_fun_x, cov_par = ego_cov_par, cov_fun = ego_cov_fun, delta = delta / 1000)
  K22 <- make_cov_matC(x = obj_fun_x, x_pred = matrix(), cov_fun = ego_cov_fun, cov_par = ego_cov_par, delta = delta / 1000)

  pred_obj_fun <- rep(obj_fun_vals[1], times = nrow(K12)) + as.numeric(K12 %*% solve(a = K22, b = obj_fun_vals - rep(obj_fun_vals[1], times = ncol(K12)))) ## objective function predictions
  obj_fun_vars <- numeric() ## objective function marginal variances
  std_vals <- numeric() ## standardized predictions
  EI <- numeric() ## expected improvement
  for(i in 1:nrow(K12))
  {
    obj_fun_vars[i] <- ego_cov_par$sigma^2 + ego_cov_par$tau^2 - K12[i,] %*% solve(a = K22, b = t(t(K12[i,])))

    if(obj_fun_vars[i] > ego_cov_par$tau^2)
    {
      std_vals[i] <- (pred_obj_fun[i] - max(obj_fun_vals)) / sqrt(obj_fun_vars[i])
      EI[i] <- sqrt(obj_fun_vars[i]) * ( std_vals[i] * pnorm(q = std_vals[i], mean = 0, sd = 1) + dnorm(x = std_vals[i], mean = 0, sd = 1) )
    }
    else{EI[i] <- 0}

  }
  iter <- TTmin
  while(iter < TTmax)
  {

    pseudo_prop <- matrix(laplace_opt$xy[which.max(EI),], nrow = 1, ncol = ncol(laplace_opt$xy))
    # print(pseudo_prop)
    # print(obj_fun_x)
    # if(
    #   !is.na(prodlim::row.match(x = as.data.frame(pseudo_prop),
    #                       table = as.data.frame(obj_fun_x)))
    #   )
    # {
    #   break
    # }

    ## meta model x values
    pseudo_xu <- rbind(xu, pseudo_prop)

    ## meta model y values
    nr_opt <- try(newtrap_sparseGP(start_vals = laplace_opt$fmax,
                               # obj_fun = obj_fun,
                               # grad_loglik_fn = grad_loglik_fn,
                               # dlog_py_dff = dlog_py_dff,
                               # d2log_py_dff = d2log_py_dff,
                               # maxit = maxit_nr,
                               # tol = tol_nr,
                               cov_par = laplace_opt$cov_par,
                               # cov_fun = cov_fun,
                               xy = xy,
                               xu = pseudo_xu,
                               # y = y,
                               mu = laplace_opt$mu,
                               muu = c(laplace_opt$muu, laplace_opt$muu[1]), delta = delta, ...))

    if(class(nr_opt) == "try-error")
    {
      while(class(nr_opt) == "try-error")
      {
        pseudo_prop <- laplace_opt$xy[sample.int(n = nrow(laplace_opt$xy), size = 1, replace = FALSE),] +
          rnorm(n = length(pseudo_prop), mean = 0, sd = laplace_opt$cov_par$l / 100)
        pseudo_xu <- rbind(xu, pseudo_prop)
        # print(pseudo_xu)

        ## meta model y values
        nr_opt <- try(newtrap_sparseGP(start_vals = laplace_opt$fmax,
                                       # obj_fun = obj_fun,
                                       # grad_loglik_fn = grad_loglik_fn,
                                       # dlog_py_dff = dlog_py_dff,
                                       # d2log_py_dff = d2log_py_dff,
                                       # maxit = maxit_nr,
                                       # tol = tol_nr,
                                       cov_par = laplace_opt$cov_par,
                                       # cov_fun = cov_fun,
                                       xy = xy,
                                       xu = pseudo_xu,
                                       # y = y,
                                       mu = laplace_opt$mu,
                                       muu = c(laplace_opt$muu, laplace_opt$muu[1]), delta = delta, ...))
      }
    }

    obj_fun_vals <- c(obj_fun_vals, nr_opt$objective_function_values[length(nr_opt$objective_function_values)])
    obj_fun_x <- rbind(obj_fun_x, pseudo_prop)

    ## optimize meta model covariance parameters
    temp_opt <- try(optim(par = log(as.numeric(ego_cov_par))[1:2], fn = meta_mod_obj_fun, gr = meta_mod_grad, method = "BFGS",
                          obj_fun_vals = obj_fun_vals, x = obj_fun_x, cov_fun = ego_cov_fun,
                          mu = obj_fun_vals[1], dcov_fun_dtheta = ego_dcov_fun_dtheta,
                          transform = TRUE, cov_par_names = names(ego_cov_par)[1:2], tau = ego_cov_par$tau, delta = delta))

    ## If optim hits values of sigma that cause infinite variance or zero length scale, then set parameter values
    ##  to something semi reasonable
    if(class(temp_opt) == "try-error")
    {
      print("Error: Optimizing meta GP was unsuccessful - likely due to extreme parameter values. Setting covariance parameters manually.")
      temp_opt <- list()
      temp_opt$par <- c(log(max(obj_fun_vals) - min(obj_fun_vals)), log( (max(obj_fun_x) - min(obj_fun_x))/1000 ))
    }
    cov_par_names <- names(ego_cov_par)
    ego_cov_par <- as.list(c(exp(temp_opt$par), ego_cov_par$tau))
    names(ego_cov_par) <- cov_par_names

    K12 <- make_cov_matC(x = laplace_opt$xy, x_pred = obj_fun_x, cov_par = ego_cov_par, cov_fun = ego_cov_fun, delta = delta / 1000)
    K22 <- make_cov_matC(x = obj_fun_x, x_pred = matrix(), cov_fun = ego_cov_fun, cov_par = ego_cov_par, delta = delta / 1000)

    pred_obj_fun <- rep(obj_fun_vals[1], times = nrow(K12)) + as.numeric(K12 %*% solve(a = K22, b = obj_fun_vals - rep(obj_fun_vals[1], times = ncol(K12)))) ## objective function predictions
    obj_fun_vars <- numeric() ## objective function marginal variances
    std_vals <- numeric() ## standardized predictions
    EI <- numeric() ## expected improvement
    for(i in 1:nrow(K12))
    {
      obj_fun_vars[i] <- ego_cov_par$sigma^2 + ego_cov_par$tau^2 - K12[i,] %*% solve(a = K22, b = t(t(K12[i,])))

      if(obj_fun_vars[i] > ego_cov_par$tau^2)
      {
        std_vals[i] <- (pred_obj_fun[i] - max(obj_fun_vals)) / sqrt(obj_fun_vars[i])
        EI[i] <- sqrt(obj_fun_vars[i]) * ( std_vals[i] * pnorm(q = std_vals[i], mean = 0, sd = 1) + dnorm(x = std_vals[i], mean = 0, sd = 1) )
      }
      else{EI[i] <- 0}

    }
    # for(i in 1:nrow(K12))
    # {
    #   obj_fun_vars[i] <- ego_cov_par$sigma^2 + ego_cov_par$tau^2 + delta - K12[i,] %*% solve(a = K22, b = t(t(K12[i,])))
    # }

    # std_vals <- (pred_obj_fun - max(obj_fun_vals)) / sqrt(obj_fun_vars)
    # EI <- sqrt(obj_fun_vars) * ( std_vals * pnorm(q = std_vals, mean = 0, sd = 1) + dnorm(x = std_vals, mean = 0, sd = 1) )

    iter <- iter + 1
  }

  # return(matrix(nrow = 1, ncol = ncol(xu),
  #               data = matrix(obj_fun_x[(nrow(xu) + 1):nrow(obj_fun_x),], ncol = ncol(obj_fun_x))[which.max(obj_fun_vals[(nrow(xu) + 1):nrow(obj_fun_x)]), ]))
  return(matrix(nrow = 1, ncol = ncol(xu),
                data = obj_fun_x[which.max(obj_fun_vals),]))

}


## knot proposal based on EGO for Gaussian data
## knot prop EGO
#'@export
knot_prop_ego_norm <- function(norm_opt,
                          # obj_fun,
                          # grad_loglik_fn,
                          # d2log_py_dff,
                          # maxit_nr = 2000,
                          # tol_nr = 1e-5,
                          # y,
                          opt = list(), ...)
{
  ## laplace_opt is a list of output returned from laplace_grad_ascent
  ## obj_fun is the objective function
  ## ... should contain the covariance function to be used to make proposals
  ##      call this function ego_cov_fun(x1, x2, cov_par)
  ## ... Must also contain the argument TTmin denoting the minimum number of allowed objective function evaluations
  ## ... Must also contain the argument TTmin denoting the maximum number of allowed objective function evaluations
  ## ... Must also contain the argument ego_cov_par: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_cov_fun: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_dcov_fun_dtheta: a list of functions which corresponds
  ##        (in order) to the covariance parameters in
  ##        cov_par. These functions compute the derivatives of the
  ##        covariance function wrt the particular covariance parameter.
  ## ... Must also contain the argument predict_laplace: to get predictive variances for the first prediction

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 0.8,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = Inf, "maxit_nr" = 1000, "delta" = 1e-6,
                    "tol_knot" = 1e-3, "maxknot" = 30, "tol_nr" = 1e-5, "TTmin" = 10, "TTmax" = 40,
                    "ego_cov_par" = list("sigma" = 3, "l" = 1, "tau" = 0), "ego_cov_fun" = "exp",
                    "ego_cov_fun_dtheta" = list("sigma" = dexp_dsigma,
                                                "l" = dexp_dl))

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: invalid or unnecessary argument to opt(). Proceeding with fingers crossed.")
        next
      }

      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  optim_method = opt_master$optim_method
  optim_par = list("decay" = opt_master$decay,
                   "epsilon" = opt_master$epsilon,
                   "eta" = opt_master$eta,
                   "learn_rate" = opt_master$learn_rate)
  maxit = opt_master$maxit
  obj_tol = opt_master$obj_tol
  grad_tol = opt_master$grad_tol
  delta <- opt_master$delta
  tol_knot <- opt_master$tol_knot
  maxknot <- opt_master$maxknot
  TTmin <- opt_master$TTmin
  TTmax <- opt_master$TTmax
  ego_cov_par <- opt_master$ego_cov_par
  ego_cov_fun <- opt_master$ego_cov_fun
  ego_dcov_fun_dtheta <- opt_master$ego_cov_fun_dtheta


  args <- list(...)
  y <- args$y
  obj_fun <- args$obj_fun
  predict_laplace <- args$predict_laplace
  cov_fun <- args$cov_fun

  ## store current knot locations
  xu <- norm_opt$xu
  cov_par <- norm_opt$cov_par
  xy <- norm_opt$xy

  ## store objective function values
  # obj_fun_vals <- rep(norm_opt$obj_fun[length(norm_opt$obj_fun)], times = nrow(xu))
  obj_fun_vals <- numeric()

  ## first set of proposals for the new knot
  # pred_vars <- as.numeric(diag(predict_laplace(u_mean = norm_opt$u_mean,
  #                          u_var = norm_opt$u_var,
  #                          xu = norm_opt$xu,
  #                          x_pred = norm_opt$xy,
  #                          cov_fun = norm_opt$cov_fun,
  #                          cov_par = norm_opt$cov_par,
  #                          mu = norm_opt$mu,
  #                          muu = norm_opt$muu)$pred_var))


  pseudo_prop <- matrix(norm_opt$xy[sample.int(n = nrow(norm_opt$xy), size = TTmin, replace = FALSE),],
                        ncol = ncol(norm_opt$xy), nrow = TTmin)

  ## meta model x values
  for(i in 1:nrow(pseudo_prop))
  {
    pseudo_xu <- rbind(xu, pseudo_prop[i,])

    ## Create data GP covariance matrices
    if(norm_opt$cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(norm_opt$xy), sep = "")
      Sigma12 <- make_cov_mat_ardC(x = norm_opt$xy, x_pred = pseudo_xu,
                               cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun,
                               delta = delta, lnames = lnames)
      Sigma22 <- make_cov_mat_ardC(x = pseudo_xu, x_pred = matrix(),
                               cov_par = norm_opt$cov_par,
                               cov_fun = norm_opt$cov_fun,
                               delta = delta,
                               lnames = lnames) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }
    else{
      Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                               cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
      Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                               cov_par = norm_opt$cov_par,
                               cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }
    ## create Z
    Z2 <- solve(a = Sigma22, b = t(Sigma12))
    Z3 <- Sigma12 * t(Z2)
    Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
    Z <- norm_opt$cov_par$sigma^2 + norm_opt$cov_par$tau^2 + delta - Z4

    ## meta model y values
    obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y))
    if(class(obj_fun_eval) == "try-error")
    {
      while(class(obj_fun_eval) == "try-error")
      {
        pseudo_prop[i,] <- norm_opt$xy[sample.int(n = nrow(norm_opt$xy), size = 1, replace = FALSE),] +
          rnorm(n = length(pseudo_prop[i,]), mean = 0, sd = norm_opt$cov_par$l / 100)
        pseudo_xu <- rbind(xu, pseudo_prop[i,])

        ## Create data GP covariance matrices
        if(norm_opt$cov_fun == "ard")
        {
          lnames <- paste("l", 1:ncol(norm_opt$xy), sep = "")
          Sigma12 <- make_cov_mat_ardC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun,
                                   delta = delta, lnames = lnames)
          Sigma22 <- make_cov_mat_ardC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta,
                                   lnames = lnames) -
            as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }
        else{
          Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
          Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }

        ## create Z
        Z2 <- solve(a = Sigma22, b = t(Sigma12))
        Z3 <- Sigma12 * t(Z2)
        Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
        Z <- norm_opt$cov_par$sigma^2 + norm_opt$cov_par$tau^2 + delta - Z4

        ## meta model y values
        obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y))
      }
    }

    obj_fun_vals <- c(obj_fun_vals, obj_fun_eval)
  }
  # obj_fun_x <- rbind(xu, pseudo_prop)
  obj_fun_x <- rbind(pseudo_prop)


  ## get meta model parameters

  ## meta model GP negative log likelihood function
  meta_mod_obj_fun <- function(cov_par, ...)
  {
    # y, x, cov_fun, mu
    args <- list(...)
    x <- args$x
    obj_fun_vals <- args$obj_fun_vals
    cov_fun <- args$cov_fun
    mu <- args$mu
    cov_par_names <- args$cov_par_names
    tau <- args$tau
    delta <- args$delta

    cov_par <- as.list(exp(cov_par))
    names(cov_par) <- cov_par_names
    cov_par$tau <- tau

    # print(cov_par)
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x), sep = "")
      Sig <- make_cov_mat_ardC(x = x, x_pred = matrix(),
                               cov_fun = cov_fun,
                               cov_par = cov_par,
                               delta = delta, lnames = lnames)

    }
    else{
      Sig <- make_cov_matC(x = x, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)

    }
    # print(Sig)

    L <- t(chol(Sig))
    mu <- rep(mu,times = length(obj_fun_vals))

    temp <- -1 * (-1 * sum(log(diag(L))) - 1/2 * t(solve(a = L, b = obj_fun_vals - mu)) %*% solve(a = L, b = obj_fun_vals - mu))

    return(temp)
  }

  ## function to create d(Sigma)/dtheta
  dsig_dtheta <- function(cov_par, par_name, dcov_fun_dtheta, x, transform = TRUE)
  {
    mat <- matrix(nrow = nrow(x), ncol = nrow(x))
    for(i in 1:nrow(mat))
    {
      for(j in 1:ncol(mat))
      {
        temp <- eval(
          substitute(expr = dcov_fun_dtheta$par(x1 = a, x2 = b, cov_par = c, transform = d),
                     env = list("a" = x[i,], "b" = x[j,], "c" = cov_par, "par" = par_name, "d" = transform))
        )
        mat[i,j] <- temp$derivative
      }
    }
    return(mat)
  }

  ## meta model gradient function
  meta_mod_grad <- function(cov_par, ...)
  {
    # y, x, cov_fun, mu, dcov_fun_dtheta, transform = TRUE
    args <- list(...)
    x <- args$x
    obj_fun_vals <- args$obj_fun_vals
    cov_fun <- args$cov_fun
    dcov_fun_dtheta <- args$dcov_fun_dtheta
    transform <- args$transform
    cov_par <- as.list(exp(cov_par))
    cov_par_names <- args$cov_par_names
    names(cov_par) <- cov_par_names
    mu <- args$mu
    mu <- rep(mu, times = nrow(x))
    tau <- args$tau
    cov_par$tau <- tau
    delta <- args$delta




    ## store gradient values
    grad <- numeric(length(dcov_fun_dtheta))

    ## create current covariancne matrix
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x), sep = "")
      Sig <- make_cov_mat_ardC(x = x, x_pred = matrix(),
                               cov_fun = cov_fun,
                               cov_par = cov_par,
                               delta = delta, lnames = lnames)

    }
    else{
      Sig <- make_cov_matC(x = x, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)

    }
    L <- t(chol(Sig))

    ## loop through each covariance parameter
    for(p in 1:length(dcov_fun_dtheta))
    {
      par_name <- names(cov_par)[p]

      ## next compute dSigma/dtheta_j
      dSigma_dtheta <- dsig_dtheta(par_name = par_name,
                                   cov_par = cov_par,
                                   dcov_fun_dtheta = dcov_fun_dtheta,
                                   x = x,
                                   transform = transform)
      grad[p] <- -(
        -1/2 * sum(diag( solve(a = t(L), b = solve(a = L, b = dSigma_dtheta)) )) +
          1/2 * t(solve(a = L, b = obj_fun_vals - mu)) %*% solve(a = L, b = dSigma_dtheta) %*% solve(a = t(L), b = solve(a = L, b = obj_fun_vals - mu))
      )
    }

    return(grad)
  }

  ## optimize meta model covariance parameters
  temp_opt <- try(optim(par = log(as.numeric(ego_cov_par))[1:2], fn = meta_mod_obj_fun, gr = meta_mod_grad, method = "BFGS",
                        obj_fun_vals = obj_fun_vals, x = obj_fun_x, cov_fun = ego_cov_fun,
                        mu = obj_fun_vals[1], dcov_fun_dtheta = ego_dcov_fun_dtheta,
                        transform = TRUE, cov_par_names = names(ego_cov_par)[1:2], tau = ego_cov_par$tau, delta = delta))

  ## If optim hits values of sigma that cause infinite variance or zero length scale, then set parameter values
  ##  to something semi reasonable
  if(class(temp_opt) == "try-error")
  {
    print("Error: Optimizing meta GP was unsuccessful - likely due to extreme parameter values. Setting covariance parameters manually.")
    temp_opt <- list()
    temp_opt$par <- c(log(max(obj_fun_vals) - min(obj_fun_vals)), log( (max(obj_fun_x) - min(obj_fun_x))/1000 ))
  }

  cov_par_names <- names(ego_cov_par)
  ego_cov_par <- as.list(c(exp(temp_opt$par), ego_cov_par$tau))
  names(ego_cov_par) <- cov_par_names


  K12 <- make_cov_matC(x = norm_opt$xy, x_pred = obj_fun_x, cov_par = ego_cov_par, cov_fun = ego_cov_fun, delta = delta / 1000)
  K22 <- make_cov_matC(x = obj_fun_x, x_pred = matrix(), cov_fun = ego_cov_fun, cov_par = ego_cov_par, delta = delta / 1000)

  pred_obj_fun <- rep(obj_fun_vals[1], times = nrow(K12)) + as.numeric(K12 %*% solve(a = K22, b = obj_fun_vals - rep(obj_fun_vals[1], times = ncol(K12)))) ## objective function predictions
  obj_fun_vars <- numeric() ## objective function marginal variances
  std_vals <- numeric() ## standardized predictions
  EI <- numeric() ## expected improvement
  for(i in 1:nrow(K12))
  {
    obj_fun_vars[i] <- ego_cov_par$sigma^2 + ego_cov_par$tau^2 - K12[i,] %*% solve(a = K22, b = t(t(K12[i,])))

    if(obj_fun_vars[i] > ego_cov_par$tau^2)
    {
      std_vals[i] <- (pred_obj_fun[i] - max(obj_fun_vals)) / sqrt(obj_fun_vars[i])
      EI[i] <- sqrt(obj_fun_vars[i]) * ( std_vals[i] * pnorm(q = std_vals[i], mean = 0, sd = 1) + dnorm(x = std_vals[i], mean = 0, sd = 1) )
    }
    else{EI[i] <- 0}

  }
  iter <- TTmin
  while(iter < TTmax)
  {

    pseudo_prop <- matrix(norm_opt$xy[which.max(EI),], nrow = 1, ncol = ncol(norm_opt$xy))
    # print(pseudo_prop)
    # print(obj_fun_x)
    # if(
    #   !is.na(prodlim::row.match(x = as.data.frame(pseudo_prop),
    #                             table = as.data.frame(obj_fun_x)))
    # )
    # {
    #   break
    # }

    ## meta model x values
    pseudo_xu <- rbind(xu, pseudo_prop)

    ## Create data GP covariance matrices
    if(norm_opt$cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(norm_opt$xy), sep = "")
      Sigma12 <- make_cov_mat_ardC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun,
                                   delta = delta, lnames = lnames)
      Sigma22 <- make_cov_mat_ardC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta,
                                   lnames = lnames) -
        as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }
    else{
      Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                               cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
      Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                               cov_par = norm_opt$cov_par,
                               cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }
    ## create Z
    Z2 <- solve(a = Sigma22, b = t(Sigma12))
    Z3 <- Sigma12 * t(Z2)
    Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
    Z <- norm_opt$cov_par$sigma^2 + norm_opt$cov_par$tau^2 + delta - Z4

    ## meta model y values
    obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y))

    if(class(obj_fun_eval) == "try-error")
    {
      while(class(obj_fun_eval) == "try-error")
      {
        pseudo_prop <- norm_opt$xy[sample.int(n = nrow(norm_opt$xy), size = 1, replace = FALSE),] +
          rnorm(n = length(pseudo_prop), mean = 0, sd = norm_opt$cov_par$l / 100)
        pseudo_xu <- rbind(xu, pseudo_prop)

        ## Create data GP covariance matrices
        if(norm_opt$cov_fun == "ard")
        {
          lnames <- paste("l", 1:ncol(norm_opt$xy), sep = "")
          Sigma12 <- make_cov_mat_ardC(x = norm_opt$xy, x_pred = pseudo_xu,
                                       cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun,
                                       delta = delta, lnames = lnames)
          Sigma22 <- make_cov_mat_ardC(x = pseudo_xu, x_pred = matrix(),
                                       cov_par = norm_opt$cov_par,
                                       cov_fun = norm_opt$cov_fun, delta = delta,
                                       lnames = lnames) -
            as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }
        else{
          Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
          Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }
        ## create Z
        Z2 <- solve(a = Sigma22, b = t(Sigma12))
        Z3 <- Sigma12 * t(Z2)
        Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
        Z <- norm_opt$cov_par$sigma^2 + norm_opt$cov_par$tau^2 + delta - Z4

        ## meta model y values
        obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y))
      }
    }

    obj_fun_vals <- c(obj_fun_vals, obj_fun_eval)

    # obj_fun_vals <- c(obj_fun_vals, norm_opt$objective_function_values[length(norm_opt$objective_function_values)])
    obj_fun_x <- rbind(obj_fun_x, pseudo_prop)

    ## optimize meta model covariance parameters
    temp_opt <- try(optim(par = log(as.numeric(ego_cov_par))[1:2], fn = meta_mod_obj_fun, gr = meta_mod_grad, method = "BFGS",
                          obj_fun_vals = obj_fun_vals, x = obj_fun_x, cov_fun = ego_cov_fun,
                          mu = obj_fun_vals[1], dcov_fun_dtheta = ego_dcov_fun_dtheta,
                          transform = TRUE, cov_par_names = names(ego_cov_par)[1:2], tau = ego_cov_par$tau, delta = delta))

    ## If optim hits values of sigma that cause infinite variance or zero length scale, then set parameter values
    ##  to something semi reasonable
    if(class(temp_opt) == "try-error")
    {
      print("Error: Optimizing meta GP was unsuccessful - likely due to extreme parameter values. Setting covariance parameters manually.")
      temp_opt <- list()
      temp_opt$par <- c(log(max(obj_fun_vals) - min(obj_fun_vals)), log( (max(obj_fun_x) - min(obj_fun_x))/1000 ))
    }
    cov_par_names <- names(ego_cov_par)
    ego_cov_par <- as.list(c(exp(temp_opt$par), ego_cov_par$tau))
    names(ego_cov_par) <- cov_par_names

    K12 <- make_cov_matC(x = norm_opt$xy, x_pred = obj_fun_x, cov_par = ego_cov_par, cov_fun = ego_cov_fun, delta = delta / 1000)
    K22 <- make_cov_matC(x = obj_fun_x, x_pred = matrix(), cov_fun = ego_cov_fun, cov_par = ego_cov_par, delta = delta / 1000)

    pred_obj_fun <- rep(obj_fun_vals[1], times = nrow(K12)) + as.numeric(K12 %*% solve(a = K22, b = obj_fun_vals - rep(obj_fun_vals[1], times = ncol(K12)))) ## objective function predictions
    obj_fun_vars <- numeric() ## objective function marginal variances
    std_vals <- numeric() ## standardized predictions
    EI <- numeric() ## expected improvement
    for(i in 1:nrow(K12))
    {
      obj_fun_vars[i] <- ego_cov_par$sigma^2 + ego_cov_par$tau^2 - K12[i,] %*% solve(a = K22, b = t(t(K12[i,])))

      if(obj_fun_vars[i] > ego_cov_par$tau^2)
      {
        std_vals[i] <- (pred_obj_fun[i] - max(obj_fun_vals)) / sqrt(obj_fun_vars[i])
        EI[i] <- sqrt(obj_fun_vars[i]) * ( std_vals[i] * pnorm(q = std_vals[i], mean = 0, sd = 1) + dnorm(x = std_vals[i], mean = 0, sd = 1) )
      }
      else{EI[i] <- 0}

    }

    # std_vals <- (pred_obj_fun - max(obj_fun_vals)) / sqrt(obj_fun_vars)
    # EI <- sqrt(obj_fun_vars) * ( std_vals * pnorm(q = std_vals, mean = 0, sd = 1) + dnorm(x = std_vals, mean = 0, sd = 1) )
    iter <- iter + 1
  }

  # return(matrix(nrow = 1, ncol = ncol(xu),
  #               data = matrix(obj_fun_x[(nrow(xu) + 1):nrow(obj_fun_x),], ncol = ncol(obj_fun_x))[which.max(obj_fun_vals[(nrow(xu) + 1):nrow(obj_fun_x)]), ]))
  return(matrix(nrow = 1, ncol = ncol(xu),
                data = obj_fun_x[which.max(obj_fun_vals),]))

}

############################################################################
## knot proposal function taking the best of a random sample
############################################################################
## knot prop random (non-Gaussian data)
#'@export
knot_prop_random <- function(laplace_opt,
                             # trace_term_fun = trace_term_fun,
                          # obj_fun,
                          # grad_loglik_fn,
                          # d2log_py_dff,
                          # maxit_nr = 2000,
                          # tol_nr = 1e-5,
                          # y,
                          opt = NA, ...)
{
  ## laplace_opt is a list of output returned from laplace_grad_ascent
  ## obj_fun is the objective function (log(p(y|f) + log(p(f|theta, xu)))
  ## d2log_py_dff is a function that computes the vector of second derivatives
  ##      of log(p(y|lambda)) wrt ff
  ## grad_loglik_fun is a function that computes the gradient of
  ##      psi = log(p(y|lambda)) + log(p(ff|theta)) wrt ff
  ## ... should contain the covariance function to be used to make proposals
  ##      call this function ego_cov_fun(x1, x2, cov_par)
  ## ... Must also contain the argument TTmin denoting the minimum number of allowed objective function evaluations
  ## ... Must also contain the argument TTmin denoting the maximum number of allowed objective function evaluations
  ## ... Must also contain the argument ego_cov_par: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_cov_fun: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_dcov_fun_dtheta: a list of functions which corresponds
  ##        (in order) to the covariance parameters in
  ##        cov_par. These functions compute the derivatives of the
  ##        covariance function wrt the particular covariance parameter.
  ## ... Must also contain the argument predict_laplace: to get predictive variances for the first prediction

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 0.8,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = 1e-1, "maxit_nr" = 1000, "delta" = 1e-6,
                    "tol_knot" = 1e-1, "maxknot" = 30, "tol_nr" = 1e-4, "TTmax" = 40)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: invalid or unnecessary argument to opt(). Proceeding with fingers crossed.")
        next
      }

      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  optim_method = opt_master$optim_method
  optim_par = list("decay" = opt_master$decay,
                   "epsilon" = opt_master$epsilon,
                   "eta" = opt_master$eta,
                   "learn_rate" = opt_master$learn_rate)
  maxit = opt_master$maxit
  maxit_nr <- opt_master$maxit_nr
  obj_tol = opt_master$obj_tol
  grad_tol = opt_master$grad_tol
  delta <- opt_master$delta
  tol_knot <- opt_master$tol_knot
  maxknot <- opt_master$maxknot
  tol_nr <- opt_master$tol_nr
  TTmin <- opt_master$TTmin
  TTmax <- opt_master$TTmax
  ego_cov_par <- opt_master$ego_cov_par
  ego_cov_fun <- opt_master$ego_cov_fun
  ego_dcov_fun_dtheta <- opt_master$ego_cov_fun_dtheta


  args <- list(...)
  obj_fun <- args$obj_fun
  d2log_py_dff <- args$d2log_py_dff
  dlog_py_dff <- args$dlog_py_dff
  grad_loglik_fun <- args$grad_loglik_fun
  predict_laplace <- args$predict_laplace
  cov_fun <- args$cov_fun

  ## store current knot locations
  xu <- laplace_opt$xu
  cov_par <- laplace_opt$cov_par
  xy <- laplace_opt$xy

  ## store objective function values
  obj_fun_vals <- rep(laplace_opt$obj_fun[length(laplace_opt$obj_fun)], times = nrow(xu))

  ## first set of proposals for the new knot
  # pred_vars <- as.numeric(diag(predict_laplace(u_mean = laplace_opt$u_mean,
  #                          u_var = laplace_opt$u_var,
  #                          xu = laplace_opt$xu,
  #                          x_pred = laplace_opt$xy,
  #                          cov_fun = laplace_opt$cov_fun,
  #                          cov_par = laplace_opt$cov_par,
  #                          mu = laplace_opt$mu,
  #                          muu = laplace_opt$muu)$pred_var))


  pseudo_prop <- matrix(laplace_opt$xy[sample.int(n = nrow(laplace_opt$xy), size = TTmax, replace = FALSE),],
                        ncol = ncol(laplace_opt$xy), nrow = TTmax)

  ## meta model x values
  for(i in 1:nrow(pseudo_prop))
  {
    pseudo_xu <- rbind(xu, pseudo_prop[i,])
    # print(pseudo_xu)

    ## meta model y values
    nr_opt <- try(newtrap_sparseGP(start_vals = laplace_opt$fmax,
                                   # obj_fun = obj_fun,
                                   # grad_loglik_fn = grad_loglik_fn,
                                   # dlog_py_dff = dlog_py_dff,
                                   # d2log_py_dff = d2log_py_dff,
                                   # maxit = maxit_nr,
                                   # tol = tol_nr,
                                   cov_par = laplace_opt$cov_par,
                                   # cov_fun = cov_fun,
                                   xy = xy,
                                   xu = pseudo_xu,
                                   # y = y,
                                   mu = laplace_opt$mu,
                                   muu = c(laplace_opt$muu, laplace_opt$muu[1]), delta = delta, ...))
    if(class(nr_opt) == "try-error")
    {
      while(class(nr_opt) == "try-error")
      {
        pseudo_prop[i,] <- laplace_opt$xy[sample.int(n = nrow(laplace_opt$xy), size = 1, replace = FALSE),] +
          rnorm(n = length(pseudo_prop[i,]), mean = 0, sd = laplace_opt$cov_par$l / 100)
        pseudo_xu <- rbind(xu, pseudo_prop[i,])
        # print(pseudo_xu)

        ## meta model y values
        nr_opt <- try(newtrap_sparseGP(start_vals = laplace_opt$fmax,
                                       # obj_fun = obj_fun,
                                       # grad_loglik_fn = grad_loglik_fn,
                                       # dlog_py_dff = dlog_py_dff,
                                       # d2log_py_dff = d2log_py_dff,
                                       # maxit = maxit_nr,
                                       # tol = tol_nr,
                                       cov_par = laplace_opt$cov_par,
                                       # cov_fun = cov_fun,
                                       xy = xy,
                                       xu = pseudo_xu,
                                       # y = y,
                                       mu = laplace_opt$mu,
                                       muu = c(laplace_opt$muu, laplace_opt$muu[1]), delta = delta, ...))
      }
    }
    obj_fun_x <- rbind(xu, pseudo_prop)

    obj_fun_vals <- c(obj_fun_vals, nr_opt$objective_function_values[length(nr_opt$objective_function_values)])
  }


  # return(matrix(nrow = 1, ncol = ncol(xu),
  #               data = matrix(obj_fun_x[(nrow(xu) + 1):nrow(obj_fun_x),], ncol = ncol(obj_fun_x))[which.max(obj_fun_vals[(nrow(xu) + 1):nrow(obj_fun_x)]), ]))
  return(matrix(nrow = 1, ncol = ncol(xu),
                data = obj_fun_x[which.max(obj_fun_vals),]))

}


## knot proposal based on random subset of the data locations
#'@export
knot_prop_random_norm <- function(norm_opt,
                               # obj_fun,
                               # grad_loglik_fn,
                               # d2log_py_dff,
                               # maxit_nr = 2000,
                               # tol_nr = 1e-5,
                               # y,
                               opt = list(), ...)
{
  ## laplace_opt is a list of output returned from laplace_grad_ascent
  ## obj_fun is the objective function
  ## ... should contain the covariance function to be used to make proposals
  ##      call this function ego_cov_fun(x1, x2, cov_par)
  ## ... Must also contain the argument TTmin denoting the minimum number of allowed objective function evaluations
  ## ... Must also contain the argument TTmin denoting the maximum number of allowed objective function evaluations
  ## ... Must also contain the argument ego_cov_par: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_cov_fun: covariance parameters for the meta model GP
  ## ... Must also contain the argument ego_dcov_fun_dtheta: a list of functions which corresponds
  ##        (in order) to the covariance parameters in
  ##        cov_par. These functions compute the derivatives of the
  ##        covariance function wrt the particular covariance parameter.
  ## ... Must also contain the argument predict_laplace: to get predictive variances for the first prediction

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 0.8,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = 1e-1, "maxit_nr" = 1000, "delta" = 1e-6,
                    "tol_knot" = 1e-1, "maxknot" = 30, "TTmax" = 40)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: invalid or unnecessary argument to opt(). Proceeding with fingers crossed.")
        # return()
        next
      }

      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  optim_method = opt_master$optim_method
  optim_par = list("decay" = opt_master$decay,
                   "epsilon" = opt_master$epsilon,
                   "eta" = opt_master$eta,
                   "learn_rate" = opt_master$learn_rate)
  maxit = opt_master$maxit
  obj_tol = opt_master$obj_tol
  grad_tol = opt_master$grad_tol
  delta <- opt_master$delta
  tol_knot <- opt_master$tol_knot
  maxknot <- opt_master$maxknot
  TTmin <- opt_master$TTmin
  TTmax <- opt_master$TTmax
  ego_cov_par <- opt_master$ego_cov_par
  ego_cov_fun <- opt_master$ego_cov_fun
  ego_dcov_fun_dtheta <- opt_master$ego_cov_fun_dtheta


  args <- list(...)
  y <- args$y
  obj_fun <- args$obj_fun
  predict_laplace <- args$predict_laplace
  cov_fun <- args$cov_fun

  ## store current knot locations
  xu <- norm_opt$xu
  cov_par <- norm_opt$cov_par
  xy <- norm_opt$xy

  ## store objective function values
  obj_fun_vals <- rep(norm_opt$obj_fun[length(norm_opt$obj_fun)], times = nrow(xu))

  ## first set of proposals for the new knot
  # pred_vars <- as.numeric(diag(predict_laplace(u_mean = norm_opt$u_mean,
  #                          u_var = norm_opt$u_var,
  #                          xu = norm_opt$xu,
  #                          x_pred = norm_opt$xy,
  #                          cov_fun = norm_opt$cov_fun,
  #                          cov_par = norm_opt$cov_par,
  #                          mu = norm_opt$mu,
  #                          muu = norm_opt$muu)$pred_var))


  pseudo_prop <- matrix(norm_opt$xy[sample.int(n = nrow(norm_opt$xy), size = TTmax, replace = FALSE),],
                        ncol = ncol(norm_opt$xy), nrow = TTmax)

  ## meta model x values
  for(i in 1:nrow(pseudo_prop))
  {
    pseudo_xu <- rbind(xu, pseudo_prop[i,])

    ## Create data GP covariance matrices
    if(norm_opt$cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(norm_opt$xy), sep = "")
      Sigma12 <- make_cov_mat_ardC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun,
                                   delta = delta, lnames = lnames)
      Sigma22 <- make_cov_mat_ardC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta,
                                   lnames = lnames) -
        as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }
    else{
      Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                               cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
      Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                               cov_par = norm_opt$cov_par,
                               cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }
    ## create Z
    Z2 <- solve(a = Sigma22, b = t(Sigma12))
    Z3 <- Sigma12 * t(Z2)
    Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
    Z <- norm_opt$cov_par$sigma^2 + norm_opt$cov_par$tau^2 + delta - Z4

    ## meta model y values
    obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y))
    if(class(obj_fun_eval) == "try-error")
    {
      while(class(obj_fun_eval == "try-error"))
      {
        pseudo_prop[i,] <- norm_opt$xy[sample.int(n = nrow(norm_opt$xy), size = 1, replace = FALSE),] +
          rnorm(n = length(pseudo_prop[i,]), mean = 0, sd = norm_opt$cov_par$l / 100)
        pseudo_xu <- rbind(xu, pseudo_prop[i,])

        ## Create data GP covariance matrices
        if(norm_opt$cov_fun == "ard")
        {
          lnames <- paste("l", 1:ncol(norm_opt$xy), sep = "")
          Sigma12 <- make_cov_mat_ardC(x = norm_opt$xy, x_pred = pseudo_xu,
                                       cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun,
                                       delta = delta, lnames = lnames)
          Sigma22 <- make_cov_mat_ardC(x = pseudo_xu, x_pred = matrix(),
                                       cov_par = norm_opt$cov_par,
                                       cov_fun = norm_opt$cov_fun, delta = delta,
                                       lnames = lnames) -
            as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }
        else{
          Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
          Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }
        ## create Z
        Z2 <- solve(a = Sigma22, b = t(Sigma12))
        Z3 <- Sigma12 * t(Z2)
        Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
        Z <- norm_opt$cov_par$sigma^2 + norm_opt$cov_par$tau^2 + delta - Z4

        ## meta model y values
        obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y))
      }
    }

    obj_fun_vals <- c(obj_fun_vals, obj_fun_eval)
  }
  obj_fun_x <- rbind(xu, pseudo_prop)

  # return(matrix(nrow = 1, ncol = ncol(xu),
  #               data = matrix(obj_fun_x[(nrow(xu) + 1):nrow(obj_fun_x),], ncol = ncol(obj_fun_x))[which.max(obj_fun_vals[(nrow(xu) + 1):nrow(obj_fun_x)]), ]))
  return(matrix(nrow = 1, ncol = ncol(xu),
                data = obj_fun_x[which.max(obj_fun_vals),]))

}
