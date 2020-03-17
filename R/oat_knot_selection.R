## One at a time knot selection function
## Gradient ascent algorithm for optimizing the Laplace approximation wrt the covariance parameters (and eventually inducing points)
## depends on Matrix and glmnet for sparse matrix operations and abind::abind for keeping track of xu values at each iteration
# library(Rcpp)
# source("covariance_function_derivatives.R")
# source("derivative_functions_of_data_likelihoods.R")
# source("gp_functions.R")
# Rcpp::sourceCpp("covariance_functionsC.cpp")
# source("laplace_approx_gradient.R")
# source("newtrap_sparseGP.R")
# source("laplace_gradient_ascent.R")
# source("laplace_approx_obj_funs.R")
#'@export
oat_knot_selection <- function(cov_par_start,
                               cov_fun,
                               dcov_fun_dtheta,
                               dcov_fun_dknot,
                               xu_start,
                               proposal_fun,
                               xy,
                               y,
                               ff,
                               grad_loglik_fn,
                               dlog_py_dff,
                               d2log_py_dff,
                               d3log_py_dff,
                               mu = NA,
                               muu_start = NA,
                               transform = TRUE,
                               obj_fun,
                               opt = list(),
                               verbose = FALSE,
                               ...)
{
  ## d2log_py_dff is a function that computes the vector of second derivatives
  ##      of log(p(y|f)) wrt f
  ## dlog_py_dff is a function that computes the vector of derivatives
  ##      of log(p(y|f)) wrt f
  ## grad_loglik_fun is a function that computes the gradient of
  ##      psi = log(p(y|lambda)) + log(p(ff|theta)) wrt ff
  ## cov_par_start is a list of the starting values of the covariance function parameters
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is the covariance function
  ## xy are the observed data locations (matrix where rows are observations)
  ## xu are the unobserved knot locations (matrix where rows are knot locations)
  ## proposal_fun is a function that proposes new knot locations to be optimized by gradient ascent
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## muu is the mean of the GP at each of the knot locations
  ## fmax is the current value for the f vector that maximizes log(p(y|f)) + log(p(f|theta, xu))
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
  ## dcov_fun_dknot is a single function that gives the derivative of the covariance function
  ##    wrt to the second input location x2, which I will always use for the knot location.
  ## transform is a logical argument that says whether to optimize on scales such that paramters are unconstrained
  ## obj_fun is the objective function
  ## learn_rate is the learning rate (step size) for the algorithm
  ## maxit is the maximum number of iterations before cutting it off
  ## tol is the absolute tolerance level which dictates
  ##      how small a difference the algorithm is allowed to make before stopping
  ##      to the objective function as well as how close the gradient is allowed to be to zero
  ## ... are argument to be passed to the functions taking derivatives of log(p(y|lambda)) wrt ff
  ##        as well as to the knot proposal function

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
                    "maxit" = 1000, "obj_tol" = 1e-2, "grad_tol" = 1, "maxit_nr" = 1000, "delta" = 1e-6,
                    "tol_knot" = 1e-1, "maxknot" = 30, "tol_nr" = 1e-5, "TTmin" = 10, "TTmax" = 30,
                    "ego_cov_par" = list("sigma" = 3, "l" = 1, "tau" = 1e-3), "ego_cov_fun" = "exp",
                    "ego_cov_fun_dtheta" = list("sigma" = dexp_dsigma,
                                                "l" = dexp_dl),
                    "cov_diff" = TRUE,
                    "chooseK" = TRUE)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print(paste( "Warning: you supplied an argument to opt() not utilized by this function: ",
        #              names(opt)[i], sep = ""))
        next
      }

      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  cov_diff <- opt_master$cov_diff
  chooseK <- opt_master$chooseK
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
  # TTmin <- opt_master$TTmin
  # TTmax <- opt_master$TTmax
  # ego_cov_par <- opt_master$ego_cov_par
  # ego_cov_fun <- opt_master$ego_cov_fun
  # ego_dcov_fun_dtheta <- opt_master$ego_dcov_fun_dtheta


  ## possibly specify default GP mean values
  if(!is.numeric(mu))
  {
    mu <- rep(mean(y), times = length(y))
  }
  if(!is.numeric(muu_start))
  {
    muu_start <- rep(mean(y), times = nrow(xu))
  }


  ## store optimized objective function values
  obj_fun_vals <- numeric()

  ## store covariance parameters at every knot optimization iteration
  cov_par_vals <- as.numeric(cov_par_start)

  ## keep track of number of GA steps at each iteration
  ga_steps <- numeric()

  xu <- xu_start

  ## optimize the objective fn wrt the covariance parameters only
  opt <- laplace_grad_ascent(cov_par_start = cov_par_start,
                             cov_fun = cov_fun,
                             dcov_fun_dtheta = dcov_fun_dtheta,
                             dcov_fun_dknot = NA, xu = xu,
                             xy = xy, ff = ff, y = y, grad_loglik_fn = grad_loglik_fn,
                             dlog_py_dff = dlog_py_dff,
                             d2log_py_dff = d2log_py_dff,
                             d3log_py_dff = d3log_py_dff,
                             mu = mu, muu =  muu_start,
                             transform = transform,
                             obj_fun = obj_fun,
                             opt = opt_master,
                             verbose = verbose, ...
  )

  ## optimize initial knot locations alongside covariance parameters
  # opt <- laplace_grad_ascent(cov_par_start = cov_par_start,
  #                            cov_fun = cov_fun,
  #                            dcov_fun_dtheta = dcov_fun_dtheta,
  #                            dcov_fun_dknot = dcov_fun_dknot, xu = xu, knot_opt = 1:nrow(xu),
  #                            xy = xy, ff = ff, y = y, grad_loglik_fn = grad_loglik_fn,
  #                            dlog_py_dff = dlog_py_dff,
  #                            d2log_py_dff = d2log_py_dff,
  #                            d3log_py_dff = d3log_py_dff,
  #                            mu = mu, muu =  muu_start,
  #                            transform = transform,
  #                            obj_fun = obj_fun,
  #                            opt = opt_master,
  #                            verbose = verbose, ...
  # )

  current_obj_fun <- opt$obj_fun[length(opt$obj_fun)] ## current objective function value
  obj_fun_vals[1] <- opt$obj_fun[length(opt$obj_fun)] ## store objective function value histories
  cov_par_vals <- rbind(cov_par_vals, as.numeric(opt$cov_par)) ## store covariance parameter value histories
  cov_par <- opt$cov_par ## current covariance parameter list
  ga_steps[1] <- opt$iter
  ## store values of GP at that maximizes the objective fn
  fmax <- opt$fmax
  muu <- opt$muu
  xu <- opt$xu

  ## store posterior mean and variance of u|y,theta,xu
  upost <- list()
  upost[[1]] <- list("mean" = opt$u_mean, "var" = opt$u_var)

  numknots <- nrow(xu_start) ## keep track of the number of knots
  iter <- 1
  opt_success <- TRUE
  ## do the OAT knot optimization
  while(ifelse(test = chooseK, yes = numknots < maxknot && ifelse(test = length(obj_fun_vals) == 1, yes = TRUE,
                                                                  no = ifelse(test = current_obj_fun - obj_fun_vals[length(obj_fun_vals) - 1] > tol_knot,
                                                                              yes = TRUE, no = FALSE)),
               no = numknots < maxknot)
  )
  {
    iter <- iter + 1

    if(verbose == TRUE)
    {
      print(iter)
    }

    ## propose a new knot
    # knot_prop <- proposal_fun(opt, ...)
    knot_prop <- proposal_fun(opt, obj_fun = obj_fun,
                              grad_loglik_fn = grad_loglik_fn,
                              dlog_py_dff = dlog_py_dff,
                              d2log_py_dff = d2log_py_dff,
                              y = y, maxit = maxit_nr, tol = tol_nr,
                              cov_fun = cov_fun, opt = opt_master, ...)

    if(verbose == TRUE)
    {
      print(knot_prop)
    }

    ## optimize the covariance parameters and knot location
    opt_temp <- laplace_grad_ascent(cov_par_start = cov_par,
                                    cov_fun = cov_fun,
                                    dcov_fun_dtheta = dcov_fun_dtheta,
                                    dcov_fun_dknot = ifelse(test = cov_diff == TRUE,
                                                            yes = dcov_fun_dknot,
                                                            no = NA),
                                    xu = rbind(xu, knot_prop),
                                    knot_opt = ifelse(test = cov_diff == TRUE,
                                                      yes = (nrow(xu) + 1):(nrow(xu) + nrow(knot_prop)),
                                                      no = NA),
                                    xy = xy, ff = fmax, y = y, grad_loglik_fn = grad_loglik_fn,
                                    dlog_py_dff = dlog_py_dff,
                                    d2log_py_dff = d2log_py_dff,
                                    d3log_py_dff = d3log_py_dff,
                                    mu = mu, muu = c(muu, muu[1]),
                                    transform = transform,
                                    obj_fun = obj_fun,
                                    opt = opt_master,
                                    verbose = verbose, ...
    )

    opt_success <- FALSE
    if(class(opt_temp) != "try-error")
    {
      opt <- opt_temp
      opt_success <- TRUE
    }

    ga_steps[iter] <- opt$iter
    current_obj_fun <- opt$obj_fun[length(opt$obj_fun)]
    obj_fun_vals[iter] <- opt$obj_fun[length(opt$obj_fun)]
    if(ifelse(test = chooseK,
              yes = current_obj_fun - obj_fun_vals[length(obj_fun_vals) - 1] > tol_knot,
              no = TRUE)
    )
    {
      cov_par_vals <- rbind(cov_par_vals, as.numeric(opt$cov_par))

      ## store values of GP at that maximizes the objective fn
      fmax <- opt$fmax

      ## store posterior mean and variance of u|y,theta,xu
      upost[[iter]] <- list("mean" = opt$u_mean, "var" = opt$u_var)

      xu <- opt$xu
      muu <- c(muu, muu[1])
      cov_par <- opt$cov_par
      numknots <- nrow(xu) ## keep track of the number of knots
    }


  }

  if(numknots == maxknot)
  {
    return(list("cov_par" = cov_par,
                "cov_fun" = cov_fun,
                "xu" = xu,
                "fmax" = fmax,
                "obj_fun" = obj_fun_vals,
                "u_mean" = upost[[iter]]$mean,
                "muu" = muu,
                "mu" = mu,
                "u_var" = upost[[iter]]$var,
                "cov_par_history" = cov_par_vals,
                "ga_steps" = ga_steps,
                "upost" = upost,
                "success" = opt_success))
  }
  else{

    return(list("cov_par" = cov_par,
                "cov_fun" = cov_fun,
                "xu" = xu,
                "fmax" = fmax,
                "obj_fun" = obj_fun_vals,
                "u_mean" = upost[[iter - 1]]$mean,
                "muu" = muu,
                "mu" = mu,
                "u_var" = upost[[iter - 1]]$var,
                "cov_par_history" = cov_par_vals,
                "ga_steps" = ga_steps,
                "u_post" = upost,
                "success" = opt_success))
  }

}


## OAT knot selection for Gaussian data
#'@export
oat_knot_selection_norm <- function(cov_par_start,
                                    cov_fun,
                                    dcov_fun_dtheta,
                                    dcov_fun_dknot,
                                    obj_fun,
                                    xu_start,
                                    proposal_fun,
                                    xy,
                                    y,
                                    ff = NA,
                                    mu,
                                    muu_start,
                                    transform = TRUE,
                                    opt = list(),
                                    verbose = FALSE,
                                    ...)
{
  ## cov_par_start is a list of the starting values of the covariance function parameters
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is the covariance function
  ## xy are the observed data locations (matrix where rows are observations)
  ## xu are the unobserved knot locations (matrix where rows are knot locations)
  ## proposal_fun is a function that proposes new knot locations to be optimized by gradient ascent
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## muu is the mean of the GP at each of the knot locations
  ## fmax is the current value for the f vector that maximizes log(p(y|f)) + log(p(f|theta, xu))
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
  ## dcov_fun_dknot is a single function that gives the derivative of the covariance function
  ##    wrt to the second input location x2, which I will always use for the knot location.
  ## transform is a logical argument that says whether to optimize on scales such that paramters are unconstrained
  ## obj_fun is the objective function
  ## learn_rate is the learning rate (step size) for the algorithm
  ## maxit is the maximum number of iterations before cutting it off
  ## tol is the absolute tolerance level which dictates
  ##      how small a difference the algorithm is allowed to make before stopping
  ##      to the objective function as well as how close the gradient is allowed to be to zero
  ## ... are argument to be passed to the functions taking derivatives of log(p(y|lambda)) wrt ff
  ##        as well as to the knot proposal function

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = 1, "maxit_nr" = 1000, "delta" = 1e-6,
                    "tol_knot" = 1e-3, "maxknot" = 30, "TTmin" = 10, "TTmax" = 40,
                    "ego_cov_par" = list("sigma" = 3, "l" = 1, "tau" = 1e-3), "ego_cov_fun" = "exp",
                    "ego_cov_fun_dtheta" = list("sigma" = dexp_dsigma,
                                                "l" = dexp_dl),
                    "cov_diff" = TRUE,
                    "chooseK" = TRUE)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print(paste( "Warning: you supplied an argument to opt() not utilized by this function: ",
        #              names(opt)[i], sep = ""))
        next
      }

      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  chooseK <- opt_master$chooseK
  cov_diff <- opt_master$cov_diff
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
  # TTmin <- opt_master$TTmin
  # TTmax <- opt_master$TTmax
  # ego_cov_par <- opt_master$ego_cov_par
  # ego_cov_fun <- opt_master$ego_cov_fun
  # ego_dcov_fun_dtheta <- opt_master$ego_dcov_fun_dtheta


  ## possibly specify default GP mean values
  if(!is.numeric(mu))
  {
    mu <- rep(mean(y), times = length(y))
  }
  if(!is.numeric(muu_start))
  {
    muu_start <- rep(mean(y), times = nrow(xu))
  }


  ## store optimized objective function values
  obj_fun_vals <- numeric()

  ## store covariance parameters at every knot optimization iteration
  cov_par_vals <- as.numeric(cov_par_start)

  ## keep track of number of GA steps at each iteration
  ga_steps <- numeric()

  xu <- xu_start
  muu <- muu_start

  ## optimize the objective fn wrt the covariance parameters only
  opt <- norm_grad_ascent(cov_par_start = cov_par_start,
                          cov_fun = cov_fun,
                          dcov_fun_dtheta = dcov_fun_dtheta,
                          dcov_fun_dknot = NA, xu = xu,
                          xy = xy, ff = NA, y = y,
                          mu = mu, muu =  muu,
                          transform = transform,
                          obj_fun = obj_fun,
                          opt = opt_master,
                          verbose = verbose, ...
  )

  current_obj_fun <- opt$obj_fun[length(opt$obj_fun)] ## current objective function value
  obj_fun_vals[1] <- opt$obj_fun[length(opt$obj_fun)] ## store objective function value histories
  cov_par_vals <- rbind(cov_par_vals, as.numeric(opt$cov_par)) ## store covariance parameter value histories
  cov_par <- opt$cov_par ## current covariance parameter list
  ga_steps[1] <- opt$iter

  ## store posterior mean and variance of u|y,theta,xu
  upost <- list()
  upost[[1]] <- list("mean" = opt$u_mean, "var" = opt$u_var)

  numknots <- nrow(xu_start) ## keep track of the number of knots
  iter <- 1
  ## do the OAT knot optimization
  while(ifelse(test = chooseK, yes = numknots < maxknot && ifelse(test = length(obj_fun_vals) == 1, yes = TRUE,
                                                                  no = ifelse(test = current_obj_fun - obj_fun_vals[length(obj_fun_vals) - 1] > tol_knot,
                                                                              yes = TRUE, no = FALSE)),
               no = numknots < maxknot)
  )
  {
    iter <- iter + 1

    if(verbose == TRUE)
    {
      print(iter)
    }

    ## propose a new knot
    # knot_prop <- proposal_fun(opt, ...)
    knot_prop <- proposal_fun(norm_opt = opt, obj_fun = obj_fun,
                              y = y, cov_fun = cov_fun, opt = opt_master, ...)

    if(verbose == TRUE)
    {
      print(knot_prop)
    }

    ## optimize the covariance parameters and knot location
    opt <- norm_grad_ascent(cov_par_start = cov_par,
                            cov_fun = cov_fun,
                            dcov_fun_dtheta = dcov_fun_dtheta,
                            dcov_fun_dknot = ifelse(test = cov_diff == TRUE,
                                                    yes = dcov_fun_dknot,
                                                    no = NA),
                            xu = rbind(xu, knot_prop),
                            knot_opt = ifelse(test = cov_diff == TRUE,
                                              yes = (nrow(xu) + 1):(nrow(xu) + nrow(knot_prop)),
                                              no = NA),
                            xy = xy, ff = NA, y = y,
                            mu = mu, muu = c(muu, muu[1]),
                            transform = transform,
                            obj_fun = obj_fun,
                            opt = opt_master,
                            verbose = verbose, ...
    )
    ga_steps[iter] <- opt$iter
    current_obj_fun <- opt$obj_fun[length(opt$obj_fun)]
    obj_fun_vals[iter] <- opt$obj_fun[length(opt$obj_fun)]


    if(ifelse(test = chooseK,
              yes = current_obj_fun - obj_fun_vals[length(obj_fun_vals) - 1] > tol_knot,
              no = TRUE)
    )
    {
      cov_par_vals <- rbind(cov_par_vals, as.numeric(opt$cov_par))

      ## store posterior mean and variance of u|y,theta,xu
      upost[[iter]] <- list("mean" = opt$u_mean, "var" = opt$u_var)

      xu <- opt$xu
      muu <- c(muu, muu[1])
      cov_par <- opt$cov_par
      numknots <- nrow(xu) ## keep track of the number of knots
    }


  }

  if(numknots == maxknot)
  {
    return(list("cov_par" = cov_par,
                "cov_fun" = cov_fun,
                "xu" = xu,
                "obj_fun" = obj_fun_vals,
                "u_mean" = upost[[iter]]$mean,
                "muu" = muu,
                "mu" = mu,
                "u_var" = upost[[iter]]$var,
                "cov_par_history" = cov_par_vals,
                "ga_steps" = ga_steps,
                "upost" = upost))
  }
  else{

    return(list("cov_par" = cov_par,
                "cov_fun" = cov_fun,
                "xu" = xu,
                "obj_fun" = obj_fun_vals,
                "u_mean" = upost[[iter - 1]]$mean,
                "muu" = muu,
                "mu" = mu,
                "u_var" = upost[[iter - 1]]$var,
                "cov_par_history" = cov_par_vals,
                "ga_steps" = ga_steps,
                "u_post" = upost))
  }

}
