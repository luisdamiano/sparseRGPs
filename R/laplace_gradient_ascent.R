# source("covariance_function_derivatives.R")
# source("derivative_functions_of_data_likelihoods.R")
# source("gp_functions.R")
# source("laplace_approx_gradient.R")
# source("newtrap_sparseGP.R")
## Gradient ascent algorithm for optimizing the Laplace approximation wrt the covariance parameters (and eventually inducing points)
## depends on Matrix and glmnet for sparse matrix operations and abind::abind for keeping track of xu values at each iteration
## depends on dlogq_dcov_par to compute gradients of Laplace Approximation
#'@export
laplace_grad_ascent <- function(cov_par_start,
                                cov_fun,
                                dcov_fun_dtheta,
                                dcov_fun_dknot,
                                knot_opt,
                                xu,
                                xy,
                                y,
                                ff,
                                grad_loglik_fn,
                                dlog_py_dff,
                                d2log_py_dff,
                                d3log_py_dff,
                                mu,
                                muu,
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
  ## cov_fun is a string specifying the covariance function ("sqexp" or "exp")
  ## xy are the observed data locations (matrix where rows are observations)
  ## knot_opt is a vector giving the knots over which to optimize, other knots are held constant
  ## xu are the unobserved knot locations (matrix where rows are knot locations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## muu is the mean of the GP at each of the knot locations
  ## ff is the starting value for the vector f
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
  ## dcov_fun_dknot is a single function that gives the derivative of the covariance function
  ##    wrt to the second input location x2, which I will always use for the knot location.
  ## transform is a logical argument that says whether to optimize on scales such that paramters are unconstrained
  ## obj_fun is the objective function
  ## learn_rate is the learning rate (step size) for the algorithm
  ## maxit is the maximum number of iterations before cutting it off
  ## delta is a small value that is added to the diagonal of covariance matrices for numerical stability
  ## tol is the absolute tolerance level which dictates
  ##      how small a difference the algorithm is allowed to make before stopping
  ##      to the objective function as well as how close the gradient is allowed to be to zero
  ## ... are argument to be passed to the functions taking derivatives of log(p(y|lambda)) wrt ff


  ## jiggle xu to make sure that knot locations aren't exactly at data locations
  # for(i in 1:nrow(xu))
  # {
  #   if(any(
  #     apply(X = xy, MARGIN = 1, FUN = function(x,y){isTRUE(all.equal(target = y, current = x))}, y = xu[i,]) == TRUE
  #   )
  #   )
  #   {
  #     xu[i,] <- xu[i,] + rnorm(n = ncol(xu), sd = 1e-6)
  #   }
  # }
  # xu <- rnorm(n = nrow(xu) * ncol(xu), sd = 1e-6)

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
             "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
             "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = Inf,
             "maxit_nr" = 1000, "delta" = 1e-6, "tol_nr" = 1e-6)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: you supplied an argument to opt() not utilized by this function.")
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
  tol_nr <- opt_master$tol_nr


  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
  }
  else{
    lnames <- character()
  }

  ## possibly specify default GP mean values
  if(!is.numeric(mu))
  {
    mu <- rep(mean(y), times = length(y))
  }
  if(!is.numeric(muu))
  {
    muu <- rep(mean(y), times = nrow(xu))
  }

  ## store objective function evaluations, gradient values, and parameter values
  ## initialize vector of objective function values
  current_cov_par <- unlist(cov_par_start) ## vector of covariance parameter values (without inducing points)
  cov_par_vals <- matrix(ncol = length(current_cov_par)) ## store the covariance parameter values
  cov_par <- as.list(current_cov_par)

  obj_fun_vals <-  numeric() ## store the objective function values
  current_obj_fun <- NA ## current objective function value

  if(is.list(dcov_fun_dtheta))
  {
    grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
    current_grad_val_theta <- NA ## current gradient value

    ## create a zero vector that will be a part of the gradient vector for the constant parameters
    # null_grad <- rep(0, length = length(current_cov_par) - length(dcov_fun_dtheta))
  }
  if(!is.list(dcov_fun_dtheta))
  {
    current_grad_val_theta <- 0 ## current gradient value
    grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
  }

  grad_knot_vals <- matrix(ncol = nrow(xu) * ncol(xu)) ## store gradient values
  current_grad_val_knot <- NA ## current gradient value

  nr_iter <- numeric() ## vector of newton raphson iterations for each GA step

  ## find the value of the GP that maximizes log(p(y|f)) + log(p(f|theta, xu))
  nr_step <- newtrap_sparseGP(start_vals = ff,
                              obj_fun = obj_fun,
                              grad_loglik_fn = grad_loglik_fn,
                              dlog_py_dff = dlog_py_dff,
                              d2log_py_dff = d2log_py_dff,
                                    maxit = maxit_nr,
                              tol = tol_nr,
                              cov_par = cov_par_start,
                              cov_fun = cov_fun,
                                    xy = xy, xu = xu,
                              y = y,
                              mu = mu,
                              muu = muu,
                              delta = delta, ...)

  nr_iter[1] <- length(nr_step$objective_function_values) ## store number of newton raphson iterations

  current_obj_fun <- nr_step$objective_function_values[length(nr_step$objective_function_values)] ## store the first objective function value log(q(y|theta, xu, fmax))
  obj_fun_vals[1] <- current_obj_fun
  cov_par_vals[1,] <- current_cov_par
  fmax <- nr_step$gp ## value of GP that maximizes the objective function
  temp_grad_eval <- dlogq_dcov_par(cov_par = cov_par_start, ## store the results of the gradient function
                              cov_fun = cov_fun,
                              dcov_fun_dtheta = dcov_fun_dtheta,
                              dcov_fun_dknot = dcov_fun_dknot,
                              knot_opt = knot_opt,
                              xu = xu,
                              xy = xy,
                              y = y,
                              ff = fmax,
                              dlog_py_dff = dlog_py_dff,
                              d2log_py_dff = d2log_py_dff,
                              d3log_py_dff = d3log_py_dff,
                              mu = mu,
                              transform = transform, delta = delta, ...)

  if(is.list(dcov_fun_dtheta))
  {
    current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
    grad_vals[1,] <- current_grad_val_theta
  }
  current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters


  ## potentially record the gradient of the knot locations
  if(is.function(dcov_fun_dknot))
  {
    current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
    grad_knot_vals[1,] <- current_grad_val_knot
    xu_vals <- xu

    ## get bounds for knots
    knot_bounds_l <- apply(X = xy, MARGIN = 2, FUN = min)
    knot_bounds_u <- apply(X = xy, MARGIN = 2, FUN = max)
    diffs <- knot_bounds_u - knot_bounds_l
    knot_bounds <- cbind(knot_bounds_l - diffs/10, knot_bounds_u + diffs/10)

    current_xu_trans <- temp_grad_eval$trans_knot
  }
  if(!is.function(dcov_fun_dknot))
  {
    current_grad_val_knot <- 0
  }


  ## run gradient ascent updates until convergence criteria is met or the maxit is reached
  iter <- 1
  # obj_fun_vals[1] <- 1e4 ## this is a line to make the objective fn check work in while loop...

  ## optimization for optim_method = "ga"
  if(optim_method == "ga")
  {
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta, current_grad_val_knot)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {

      ## update the iteration counter
      iter <- iter + 1

      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta, current_grad_val_knot))
        # print(xu)
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a GA step in transformed space
        if(is.list(dcov_fun_dtheta))
        {
          current_cov_par_trans <- current_cov_par_trans + optim_par$learn_rate * current_grad_val_theta

          ## recover the untransformed parameter values
          for(j in 1:length(cov_par))
          {
            ## transform the transformed parameters back to the original parameter space
            if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
            {
              current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
            }
            else{
              current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                  cov_par = cov_par_start,
                                                                                                                  transform = TRUE)$trans_fun(current_cov_par_trans[j])
            }

          }
          cov_par <- as.list(current_cov_par)
        }

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
          current_xu_trans <- current_xu_trans + optim_par$learn_rate * matrix(data = current_grad_val_knot,
                                                                               nrow = nrow(xu),
                                                                               ncol = ncol(xu),
                                                                               byrow = TRUE)
          xu <- apply(X = current_xu_trans, MARGIN = 1,
                      FUN = dcov_fun_dknot(x1 = 0, x2 = 0, cov_par = cov_par, transform = TRUE, bounds = knot_bounds)$trans_fun,
                      bounds = knot_bounds)

          xu <- matrix(xu, ncol = ncol(xy), byrow = TRUE)

          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }
      }
      if(transform == FALSE)
      {
        current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
        cov_par <- as.list(current_cov_par)
        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
          xu <- xu + optim_par$learn_rate* matrix(data = current_grad_val_knot,
                                              nrow = nrow(xu),
                                              ncol = ncol(xu),
                                              byrow = TRUE)
          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }
      }

      ## find the value of the GP that maximizes log(p(y|f)) + log(p(f|theta, xu))
      nr_step <- newtrap_sparseGP(start_vals = fmax,
                                  obj_fun = obj_fun,
                                  grad_loglik_fn = grad_loglik_fn,
                                  dlog_py_dff = dlog_py_dff,
                                  d2log_py_dff = d2log_py_dff,
                                  maxit = maxit_nr,
                                  tol = tol_nr,
                                  cov_par = cov_par,
                                  cov_fun = cov_fun,
                                  xy = xy, xu = xu, y = y, mu = mu, muu = muu, delta = delta, ...)

      nr_iter[iter] <- length(nr_step$objective_function_values)

      current_obj_fun <- nr_step$objective_function_values[length(nr_step$objective_function_values)] ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)
      fmax <- nr_step$gp ## value of GP that maximizes the objective function
      temp_grad_eval <- dlogq_dcov_par(cov_par = cov_par, ## store the results of the gradient function
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = knot_opt,
                                       xu = xu,
                                       xy = xy,
                                       y = y,
                                       ff = fmax,
                                       dlog_py_dff = dlog_py_dff,
                                       d2log_py_dff = d2log_py_dff,
                                       d3log_py_dff = d3log_py_dff,
                                       mu = mu,
                                       delta = delta,
                                       transform = transform, ...)

      if(is.list(dcov_fun_dtheta))
      {
        current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
        current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
        grad_vals <- rbind(grad_vals, current_grad_val_theta)
      }

      ## potentially record the gradient of the knot locations
      if(is.function(dcov_fun_dknot))
      {
        current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
        grad_knot_vals <- rbind(grad_knot_vals, current_grad_val_knot)
      }

    }
  }

  ## optimization for optim_method = "adadelta"
  sg2_theta <- 0
  sg2_knot <- 0
  sd2_theta <- 0
  sd2_knot <- 0
  if(optim_method == "adadelta")
  {
    if(is.list(dcov_fun_dtheta))
    {
      ## store the decaying sum of squared gradients
      sg2_theta <- rep(0, times = length(current_cov_par))

      ## store the decaying sum of updates
      sd2_theta <- rep(0, times = length(current_cov_par))

      sign_change_theta <-  rep(0, times = length(current_cov_par))
    }
    if(is.function(dcov_fun_dknot))
    {
      ## store the decaying sum of squared gradients
      sg2_knot <- rep(0, times = nrow(xu)*ncol(xu))

      ## store the decaying sum of updates
      sd2_knot <- rep(0, times = nrow(xu)*ncol(xu))

      sign_change_knot <- rep(0, times = nrow(xu)*ncol(xu))
    }

    if(verbose == TRUE)
    {
      print(paste("iteration ", iter, sep = ""))
      print(c(current_grad_val_theta, current_grad_val_knot))
      # print(current_cov_par)
      # print(xu)
      #
      # print(current_obj_fun)
    }

    ## do the optimization
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta, current_grad_val_knot)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {
      ## update the iteration counter
      iter <- iter + 1

      if(transform == TRUE)
      {
        ## take a adadelta step in transformed space
        if(is.list(dcov_fun_dtheta))
        {
          ## accumulate the gradient
          sg2_theta <- optim_par$decay * sg2_theta + (1 - optim_par$decay) * current_grad_val_theta^2

          ## calculate adapted step size
          # delta_theta <- (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta
          delta_theta <- ((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta

          # delta_theta <- optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))) * current_grad_val_theta

          # if(verbose == TRUE)
          # {
          #   print("adaptive step size theta")
          #   # print(optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))))
          #   print(((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))))
          #   # print(((1/optim_par$eta)^(sign_change_theta)))
          # }

          ## accumulate the step size
          sd2_theta <- optim_par$decay * sd2_theta + (1 - optim_par$decay) * delta_theta^2

          ## update the parameters
          current_cov_par_trans <- current_cov_par_trans + delta_theta
        }

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )

          ## accumulate the gradient
          sg2_knot <- optim_par$decay * sg2_knot + (1 - optim_par$decay) * current_grad_val_knot^2

          ## calculate adapted step size
          # delta_knot <- (sqrt(sd2_knot + rep(optim_par$epsilon, times = length(sd2_knot))) / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot

          delta_knot <- ((1/optim_par$eta)^(sign_change_knot)) * (sqrt(sd2_knot + rep(optim_par$epsilon, times = length(sd2_knot))) / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot
          # delta_knot <- (optim_par$learn_rate / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot

          # if(verbose == TRUE)
          # {
          #   #print("adaptive step size knot")
          #   #print(((1/optim_par$eta)^(sign_change_knot)) * (sqrt(sd2_knot + rep(optim_par$epsilon, times = length(sd2_knot))) / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))))
          #   # print(((1/optim_par$eta)^(sign_change_knot)))
          # }

          ## accumulate the step size
          sd2_knot <- optim_par$decay * sd2_knot + (1 - optim_par$decay) * delta_knot^2

          ## update knot values
          current_xu_trans <- current_xu_trans + matrix(data = delta_knot,
                                                        nrow = nrow(xu),
                                                        ncol = ncol(xu),
                                                        byrow = TRUE)
          xu <- apply(X = current_xu_trans, MARGIN = 1,
                      FUN = dcov_fun_dknot(x1 = 0, x2 = 0, cov_par = cov_par, transform = TRUE, bounds = knot_bounds)$trans_fun,
                      bounds = knot_bounds)

          xu <- matrix(xu, ncol = ncol(xy), byrow = TRUE)

          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }

        if(is.list(dcov_fun_dtheta))
        {
          ## recover the untransformed parameter values
          for(j in 1:length(cov_par))
          {
            ## transform the transformed parameters back to the original parameter space
            if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
            {
              current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
            }
            else{
              current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                  cov_par = cov_par_start,
                                                                                                                  transform = TRUE)$trans_fun(current_cov_par_trans[j])
            }

          }
          cov_par <- as.list(current_cov_par)

        }
      }
      ## This commented code does gradient ascent... not modified adadelta
      ##  I should probably cut out the transform argument and mandate a transformation to the whole real line
      # if(transform == FALSE)
      # {
      #   current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
      #   cov_par <- as.list(current_cov_par)
      #   if(is.function(dcov_fun_dknot))
      #   {
      #     ## Note that the gradient vector for the knots looks like this
      #     ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
      #     xu <- xu + optim_par$learn_rate * matrix(data = current_grad_val_knot,
      #                                         nrow = nrow(xu),
      #                                         ncol = ncol(xu),
      #                                         byrow = TRUE)
      #     xu_vals <- abind::abind(xu_vals, xu, along = 3)
      #   }
      # }

      ## find the value of the GP that maximizes log(p(y|f)) + log(p(f|theta, xu))
      nr_step <- newtrap_sparseGP(start_vals = fmax,
                                  obj_fun = obj_fun,
                                  grad_loglik_fn = grad_loglik_fn,
                                  dlog_py_dff = dlog_py_dff,
                                  d2log_py_dff = d2log_py_dff,
                                  maxit = maxit_nr,
                                  tol = tol_nr,
                                  cov_par = cov_par,
                                  cov_fun = cov_fun,
                                  xy = xy, xu = xu, y = y, mu = mu, delta = delta, muu = muu, ...)

      nr_iter[iter] <- length(nr_step$objective_function_values)

      current_obj_fun <- nr_step$objective_function_values[length(nr_step$objective_function_values)] ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)
      fmax <- nr_step$gp ## value of GP that maximizes the objective function
      temp_grad_eval <- dlogq_dcov_par(cov_par = cov_par, ## store the results of the gradient function
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = knot_opt,
                                       xu = xu,
                                       xy = xy,
                                       y = y,
                                       ff = fmax,
                                       dlog_py_dff = dlog_py_dff,
                                       d2log_py_dff = d2log_py_dff,
                                       d3log_py_dff = d3log_py_dff,
                                       mu = mu,
                                       delta = delta,
                                       transform = transform, ...)


      if(is.list(dcov_fun_dtheta))
      {
        ## get the sign change and add it to the last one
        sign_change_theta <- (optim_par$decay * sign_change_theta +
          (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))/2)
        # sign_change_theta <- sign_change_theta +
        #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient) - sign(current_grad_val_theta))

        current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
        current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
        grad_vals <- rbind(grad_vals, current_grad_val_theta)
      }

      ## potentially record the gradient of the knot locations
      if(is.function(dcov_fun_dknot))
      {
        ## get the sign change and add it to the last one
        sign_change_knot <- (optim_par$decay * sign_change_knot +
          (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$knot_gradient)) - sign(current_grad_val_knot)))

        # sign_change_knot <- sign_change_knot +
        #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$knot_gradient)) - sign(current_grad_val_knot))

        current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
        grad_knot_vals <- rbind(grad_knot_vals, current_grad_val_knot)
      }

      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta, current_grad_val_knot))
        # print(current_cov_par)
        # print(xu)
        #
        # print(current_obj_fun)
      }

    }
  }

  if(is.function(dcov_fun_dknot))
  {
    return(list("cov_par" = cov_par,
                "cov_fun" = cov_fun,
                "xu" = xu,
                "xy" = xy,
                "mu" = mu,
                "muu" = muu,
                "fmax" = fmax,
                "iter" = iter,
                "obj_fun" = obj_fun_vals,
                "fmax" = fmax,
                "u_mean" = nr_step$u_posterior_mean,
                "u_var" = nr_step$u_posterior_variance,
                "grad" = grad_vals,
                "knot_grad" = grad_knot_vals,
                "knot_history" = xu_vals,
                "cov_par_history" = cov_par_vals))
  }
  if(!is.function(dcov_fun_dknot))
  {
    return(list("cov_par" = cov_par,
                "cov_fun" = cov_fun,
                "xu" = xu,
                "xy" = xy,
                "mu" = mu,
                "muu" = muu,
                "fmax" = fmax,
                "iter" = iter,
                "obj_fun" = obj_fun_vals,
                "fmax" = fmax,
                "u_mean" = nr_step$u_posterior_mean,
                "u_var" = nr_step$u_posterior_variance,
                "grad" = grad_vals,
                "knot_grad" = 0,
                "knot_history" = xu,
                "cov_par_history" = cov_par_vals))
  }

}

## FULL GP LAPLACE GRADIENT ASCENT
## Gradient ascent algorithm for optimizing the Laplace approximation wrt the covariance parameters (and eventually inducing points)
## depends on Matrix and glmnet for sparse matrix operations and abind::abind for keeping track of xu values at each iteration
## depends on dlogq_dcov_par to compute gradients of Laplace Approximation
#'@export
laplace_grad_ascent_full <- function(cov_par_start,
                                cov_fun,
                                dcov_fun_dtheta,
                                xy,
                                y,
                                ff,
                                grad_loglik_fn,
                                dlog_py_dff,
                                d2log_py_dff,
                                d3log_py_dff,
                                mu,
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
  ## cov_fun is a string specifying the covariance function ("sqexp" or "exp")
  ## xy are the observed data locations (matrix where rows are observations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## ff is the starting value for the vector f
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
  ## dcov_fun_dknot is a single function that gives the derivative of the covariance function
  ##    wrt to the second input location x2, which I will always use for the knot location.
  ## transform is a logical argument that says whether to optimize on scales such that paramters are unconstrained
  ## obj_fun is the objective function
  ## learn_rate is the learning rate (step size) for the algorithm
  ## maxit is the maximum number of iterations before cutting it off
  ## delta is a small value that is added to the diagonal of covariance matrices for numerical stability
  ## tol is the absolute tolerance level which dictates
  ##      how small a difference the algorithm is allowed to make before stopping
  ##      to the objective function as well as how close the gradient is allowed to be to zero
  ## ... are argument to be passed to the functions taking derivatives of log(p(y|lambda)) wrt ff

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = Inf,
                    "maxit_nr" = 1000, "delta" = 1e-6, "tol_nr" = 1e-6)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: you supplied an argument to opt() not utilized by this function.")
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
  tol_nr <- opt_master$tol_nr


  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
  }
  else{
    lnames <- character()
  }

  ## possibly specify default GP mean values
  if(!is.numeric(mu))
  {
    mu <- rep(mean(y), times = length(y))
  }

  ## store objective function evaluations, gradient values, and parameter values
  ## initialize vector of objective function values
  current_cov_par <- unlist(cov_par_start) ## vector of covariance parameter values (without inducing points)
  cov_par_vals <- matrix(ncol = length(current_cov_par)) ## store the covariance parameter values
  cov_par <- as.list(current_cov_par)

  obj_fun_vals <-  numeric() ## store the objective function values
  current_obj_fun <- NA ## current objective function value

  grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
  current_grad_val_theta <- NA ## current gradient value

  nr_iter <- numeric() ## vector of newton raphson iterations for each GA step

  ## find the value of the GP that maximizes log(p(y|f)) + log(p(f|theta, xu))
  nr_step <- newtrapGP(start_vals = ff,
                              dlog_py_dff = dlog_py_dff,
                              obj_fun = obj_fun,
                              grad_loglik_fn = grad_loglik_fn,
                              d2log_py_dff = d2log_py_dff,
                              maxit = maxit_nr,
                              tol = tol_nr,
                              cov_par = cov_par_start,
                              cov_fun = cov_fun,
                              xy = xy, y = y, mu = mu, delta = delta, ...)

  nr_iter[1] <- length(nr_step$objective_function_values) ## store number of newton raphson iterations

  current_obj_fun <- nr_step$objective_function_values[length(nr_step$objective_function_values)] ## store the first objective function value log(q(y|theta, xu, fmax))
  obj_fun_vals[1] <- current_obj_fun
  cov_par_vals[1,] <- current_cov_par
  fmax <- nr_step$gp ## value of GP that maximizes the objective function
  temp_grad_eval <- dlogq_dcov_par_full(cov_par = cov_par_start, ## store the results of the gradient function
                                   cov_fun = cov_fun,
                                   dcov_fun_dtheta = dcov_fun_dtheta,
                                   xy = xy,
                                   y = y,
                                   ff = fmax,
                                   dlog_py_dff = dlog_py_dff,
                                   d2log_py_dff = d2log_py_dff,
                                   d3log_py_dff = d3log_py_dff,
                                   mu = mu,
                                   transform = transform, delta = delta, ...)

  current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
  grad_vals[1,] <- current_grad_val_theta
  current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters

  ## run gradient ascent updates until convergence criteria is met or the maxit is reached
  iter <- 1
  # obj_fun_vals[1] <- 1e4 ## this is a line to make the objective fn check work in while loop...

  ## optimization for optim_method = "ga"
  if(optim_method == "ga")
  {
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {

      ## update the iteration counter
      iter <- iter + 1

      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta))
        # print(xu)
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a GA step in transformed space
        current_cov_par_trans <- current_cov_par_trans + optim_par$learn_rate * current_grad_val_theta

        ## recover the untransformed parameter values
        for(j in 1:length(cov_par))
        {
          ## transform the transformed parameters back to the original parameter space
          if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
          {
            current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
          }
          else{
            current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                cov_par = cov_par_start,
                                                                                                                transform = TRUE)$trans_fun(current_cov_par_trans[j])
          }

        }
        cov_par <- as.list(current_cov_par)
      }
      if(transform == FALSE)
      {
        current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
        cov_par <- as.list(current_cov_par)
      }

      ## find the value of the GP that maximizes log(p(y|f)) + log(p(f|theta, xu))
      nr_step <- newtrapGP(start_vals = fmax,
                                  obj_fun = obj_fun,
                                  dlog_py_dff = dlog_py_dff,
                                  grad_loglik_fn = grad_loglik_fn,
                                  d2log_py_dff = d2log_py_dff,
                                  maxit = maxit_nr,
                                  tol = tol_nr,
                                  cov_par = cov_par,
                                  cov_fun = cov_fun,
                                  xy = xy, y = y, mu = mu, delta = delta, ...)

      nr_iter[iter] <- length(nr_step$objective_function_values)

      current_obj_fun <- nr_step$objective_function_values[length(nr_step$objective_function_values)] ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)
      fmax <- nr_step$gp ## value of GP that maximizes the objective function
      temp_grad_eval <- dlogq_dcov_par_full(cov_par = cov_par, ## store the results of the gradient function
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       xy = xy,
                                       y = y,
                                       ff = fmax,
                                       dlog_py_dff = dlog_py_dff,
                                       d2log_py_dff = d2log_py_dff,
                                       d3log_py_dff = d3log_py_dff,
                                       mu = mu,
                                       delta = delta,
                                       transform = transform, ...)

      current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
      current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
      grad_vals <- rbind(grad_vals, current_grad_val_theta)

    }
  }

  ## optimization for optim_method = "adadelta"
  sg2_theta <- 0
  sg2_knot <- 0
  sd2_theta <- 0
  sd2_knot <- 0
  if(optim_method == "adadelta")
  {
    ## store the decaying sum of squared gradients
    sg2_theta <- rep(0, times = length(current_cov_par))

    ## store the decaying sum of updates
    sd2_theta <- rep(0, times = length(current_cov_par))

    sign_change_theta <-  rep(0, times = length(current_cov_par))

    ## do the optimization
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {
      ## update the iteration counter
      iter <- iter + 1
      #
      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta))
        # print(current_cov_par)
        # print(xu)
        #
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a adadelta step in transformed space
        ## accumulate the gradient
        sg2_theta <- optim_par$decay * sg2_theta + (1 - optim_par$decay) * current_grad_val_theta^2

        ## calculate adapted step size
        # delta_theta <- (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta
        delta_theta <- ((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta

        # delta_theta <- optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))) * current_grad_val_theta

        # if(verbose == TRUE)
        # {
        #   print("adaptive step size theta")
        #   # print(optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))))
        #   print(((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))))
        #   # print(((1/optim_par$eta)^(sign_change_theta)))
        # }

        ## accumulate the step size
        sd2_theta <- optim_par$decay * sd2_theta + (1 - optim_par$decay) * delta_theta^2

        ## update the parameters
        current_cov_par_trans <- current_cov_par_trans + delta_theta

        ## recover the untransformed parameter values
        for(j in 1:length(cov_par))
        {
          ## transform the transformed parameters back to the original parameter space
          if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
          {
            current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
          }
          else{
            current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                cov_par = cov_par_start,
                                                                                                                transform = TRUE)$trans_fun(current_cov_par_trans[j])
          }

        }
        cov_par <- as.list(current_cov_par)
      }
      ## This commented code does gradient ascent... not modified adadelta
      ##  I should probably cut out the transform argument and mandate a transformation to the whole real line
      # if(transform == FALSE)
      # {
      #   current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
      #   cov_par <- as.list(current_cov_par)
      #   if(is.function(dcov_fun_dknot))
      #   {
      #     ## Note that the gradient vector for the knots looks like this
      #     ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
      #     xu <- xu + optim_par$learn_rate * matrix(data = current_grad_val_knot,
      #                                         nrow = nrow(xu),
      #                                         ncol = ncol(xu),
      #                                         byrow = TRUE)
      #     xu_vals <- abind::abind(xu_vals, xu, along = 3)
      #   }
      # }

      ## find the value of the GP that maximizes log(p(y|f)) + log(p(f|theta, xu))
      nr_step <- newtrapGP(start_vals = fmax,
                                  obj_fun = obj_fun,
                                  dlog_py_dff = dlog_py_dff,
                                  grad_loglik_fn = grad_loglik_fn,
                                  d2log_py_dff = d2log_py_dff,
                                  maxit = maxit_nr,
                                  tol = tol_nr,
                                  cov_par = cov_par,
                                  cov_fun = cov_fun,
                                  xy = xy, y = y, mu = mu, delta = delta, ...)

      nr_iter[iter] <- length(nr_step$objective_function_values)

      current_obj_fun <- nr_step$objective_function_values[length(nr_step$objective_function_values)] ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)
      fmax <- nr_step$gp ## value of GP that maximizes the objective function
      temp_grad_eval <- dlogq_dcov_par_full(cov_par = cov_par, ## store the results of the gradient function
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       xy = xy,
                                       y = y,
                                       ff = fmax,
                                       dlog_py_dff = dlog_py_dff,
                                       d2log_py_dff = d2log_py_dff,
                                       d3log_py_dff = d3log_py_dff,
                                       mu = mu,
                                       delta = delta,
                                       transform = transform, ...)


      ## get the sign change and add it to the last one
      sign_change_theta <- (optim_par$decay * sign_change_theta +
                              (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))/2)
      # sign_change_theta <- sign_change_theta +
      #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient) - sign(current_grad_val_theta))

      current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
      current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
      grad_vals <- rbind(grad_vals, current_grad_val_theta)

    }
  }

  return(list("cov_par" = cov_par,
              "cov_fun" = cov_fun,
              "xy" = xy,
              "mu" = mu,
              "cov_fun" = cov_fun,
              "fmax" = fmax,
              "iter" = iter,
              "obj_fun" = obj_fun_vals,
              "fmax" = fmax,
              "grad" = grad_vals,
              "cov_par_history" = cov_par_vals,
              "L" = nr_step$L,
              "W" = nr_step$W,
              "grad_log_py_ff" = nr_step$grad_log_py_ff))

}

## testing gradient ascent
# y <- c(rpois(n = 2, lambda = c(2,2)))
# cov_par <- list("sigma" = 1.25, "l" = 2.27, "tau" = 1e-5)
# cov_fun <- mv_cov_fun_sqrd_exp
# dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl)
# xy <- matrix(rep(c(1,2), times = 1), nrow = 2)
# xu <- matrix(rep(c(1), times = 1), nrow = 1)
# dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl)
# ff <- log(c(2,2))
# mu <- c(0,0)
# m <- c(1,1)
#
#
## Test newton-raphson algorithm on simulated sgcp
# dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl)
# source("simulate_sgcp.R")
# set.seed(1308)
# test_sgcp <- simulate_sgcp(lambda_star = 50, gp_par = list("sigma" = 50, "l" = 0.5, "tau" = sqrt(1e-5)), dbounds = c(0,10), cov_fun = cov_fun_sqrd_exp)
# test_sgcp
# length(test_sgcp)
#
# plot(x = test_sgcp, y = rep(0, times = length(test_sgcp)))
# hist(test_sgcp, breaks = 100)
#
# cov_par <- list("sigma" = 10, "l" = 0.1, "tau" = sqrt(1e-5))
# xy <- matrix(seq(from = 0, to = 10 - 0.01, by = 0.02) + 0.01, ncol = 1)
# xu <- matrix(1:10, ncol = 1)
# y <- numeric(nrow(xy))
# for(i in 1:nrow(xy))
# {
#   y[i] <- sum(test_sgcp > (xy[i] - 0.1) & test_sgcp < (xy[i] + 0.1))
# }
# y
# m <- rep(0.02, times = length(xy))
# mu <- rep(log(sum(y)/10), times = length(xy))
# ff <- rep(log(sum(y)/10), times = length(xy))
# #
# # ## this examples takes 45 minutes
# cov_par <- list("sigma" = 4.1052, "l" = 0.70248, "tau" = sqrt(1e-5))
# system.time(test_not_oat <- laplace_grad_ascent(cov_par_start = cov_par,
#                                                 cov_fun = "sqexp",
#                                                 dcov_fun_dtheta = dcov_fun_dtheta,
#                                                 dcov_fun_dknot = NA,
#                                                 xu = xu,
#                                                 xy = xy,
#                                                 y = y,
#                                                 ff = fstart,
#                                                 grad_loglik_fn = grad_loglik_fn_pois,
#                                                 dlog_py_dff = dlog_py_dff_pois,
#                                                 d2log_py_dff = d2log_py_dff_pois,
#                                                 mu = mu,
#                                                 muu = muu,
#                                                 transform = TRUE,
#                                                 obj_fun = obj_fun_pois,
#                                                 optim_method = "adadelta",
#                                                 optim_par = list("decay" = 0.95, "epsilon" = 1e-6, "eta" = 1e4),
#                                                 maxit = 3000,
#                                                 maxit_nr = 1000,
#                                                 tol_nr = 1e-4,
#                                                 tol = 1e-3,
#                                                 verbose = TRUE,
#                                                 m = rep(area^2, times = nrow(xy))
# )
# )

# test$iter
# test$cov_par
# plot(1:(test$iter), test$grad[,1], type = "l")
# plot(1:(test$iter), test$grad[,2], type = "l")
# plot(1:(test$iter), test$obj_fun, type = "l")
# plot(1:(test$iter), test$cov_par_history[,1], type = "l")
# plot(1:(test$iter), test$cov_par_history[,2], type = "l")
#
# hist(test_sgcp, breaks = 100, ylim = c(0,max(exp(test$fmax))))
# points(x = xy, y = test$fmax, type = "l", col = "red")
# points(x = xy, y = exp(test$fmax), type = "l", col = "blue")
#
# hist(test_sgcp, breaks = 100)
# points(x = xy, y = test$fmax, type = "l", col = "red")
#
# test$cov_par
#
#
# plot(1:1000, test$grad[1:1000,1], type = "l")
# plot(1:1000, test$grad[1:1000,2])
# plot(1:1000, test$obj_fun)
# plot(1:1000, test$cov_par_history[1:1000,1])
# plot(1:1000, test$cov_par_history[1:1000,2])

## Gradient ascent algorithm for optimizing the likelihood wrt the covariance parameters (and eventually inducing points)
## depends on Matrix and glmnet for sparse matrix operations and abind::abind for keeping track of xu values at each iteration
## depends on dlogq_dcov_par to compute gradients of Laplace Approximation
#'@export
norm_grad_ascent <- function(cov_par_start,
                                cov_fun,
                                dcov_fun_dtheta,
                                dcov_fun_dknot,
                                knot_opt,
                                xu,
                                xy,
                                y,
                                mu = NA,
                                muu = NA,
                                transform = TRUE,
                                obj_fun,
                                opt = list(),
                                verbose = FALSE,
                                ...)
{
  ## cov_par_start is a list of the starting values of the covariance function parameters
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is a string specifying the covariance function ("sqexp" or "exp")
  ## xy are the observed data locations (matrix where rows are observations)
  ## knot_opt is a vector giving the knots over which to optimize, other knots are held constant
  ## xu are the unobserved knot locations (matrix where rows are knot locations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## muu is the mean of the GP at each of the knot locations
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

  ## jiggle xu to make sure that knot locations aren't exactly at data locations
  # for(i in 1:nrow(xu))
  # {
  #   if(any(
  #     apply(X = xy, MARGIN = 1, FUN = function(x,y){isTRUE(all.equal(target = y, current = x))}, y = xu[i,]) == TRUE
  #     )
  #   )
  #     {
  #     xu[i,] <- xu[i,] + rnorm(n = ncol(xu), sd = 1e-6)
  #   }
  # }

  ## extract names of optional arguments from opt
  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = Inf, "delta" = 1e-6)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: you supplied an argument to opt() not utilized by this function.")
        next
      }
      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  optim_method = opt_master$optim_method
  optim_par = list("decay" = opt_master$decay, "epsilon" = opt_master$epsilon, "eta" = opt_master$eta, "learn_rate" = opt_master$learn_rate)
  maxit = opt_master$maxit
  obj_tol = opt_master$obj_tol
  grad_tol = opt_master$grad_tol
  delta <- opt_master$delta

  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
  }
  else{
    lnames <- character()
  }

  if(!is.numeric(mu))
  {
    mu <- rep(mean(y), times = length(y))
  }
  if(!is.numeric(muu))
  {
    muu <- rep(mean(y), times = nrow(xu))
  }

  ## make sure y is numeric
  y <- as.numeric(y)

  ## store objective function evaluations, gradient values, and parameter values
  ## initialize vector of objective function values
  current_cov_par <- unlist(cov_par_start) ## vector of covariance parameter values (without inducing points)
  cov_par_vals <- matrix(ncol = length(current_cov_par)) ## store the covariance parameter values
  cov_par <- as.list(current_cov_par)

  obj_fun_vals <-  numeric() ## store the objective function values
  current_obj_fun <- NA ## current objective function value

  if(is.list(dcov_fun_dtheta))
  {
    grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
    current_grad_val_theta <- NA ## current gradient value
  }
  if(!is.list(dcov_fun_dtheta))
  {
    current_grad_val_theta <- 0 ## current gradient value
    grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
  }

  grad_knot_vals <- matrix(ncol = nrow(xu) * ncol(xu)) ## store gradient values
  current_grad_val_knot <- NA ## current gradient value

  ## current objective function value
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
    Sigma12 <- make_cov_mat_ardC(x = xy, x_pred = xu, cov_par = cov_par,
                                 cov_fun = cov_fun,
                                 delta = delta, lnames = lnames)
    Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(),
                                 cov_par = cov_par,
                                 cov_fun = cov_fun,
                                 delta = delta,
                                 lnames = lnames) -
      as.list(cov_par)$tau^2 * diag(nrow(xu))
  }
  else{
    lnames <- character()
    Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
    Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(),
                             cov_par = cov_par,
                             cov_fun = cov_fun, delta = delta) - as.list(cov_par)$tau^2 * diag(nrow(xu))
  }

  ## create Z
  Z2 <- solve(a = Sigma22, b = t(Sigma12))
  Z3 <- Sigma12 * t(Z2)
  Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
  Z <- cov_par$sigma^2 + cov_par$tau^2 + delta - Z4

  current_obj_fun <- obj_fun(mu = mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, ...) ## store the first objective function value log(q(y|theta, xu, fmax))
  obj_fun_vals[1] <- current_obj_fun
  cov_par_vals[1,] <- current_cov_par
  temp_grad_eval <- dlogp_dcov_par(cov_par = cov_par, cov_fun = cov_fun, ## store the results of gradient function
                                   dcov_fun_dtheta = dcov_fun_dtheta,
                                   dcov_fun_dknot = dcov_fun_dknot,
                                   knot_opt = knot_opt,
                                   xu = xu, xy = xy,
                                   y = y, ff = NA, mu = mu, delta = delta, transform = transform)

  if(is.list(dcov_fun_dtheta))
  {
    current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
    grad_vals[1,] <- current_grad_val_theta
  }
  current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters


  ## potentially record the gradient of the knot locations
  if(is.function(dcov_fun_dknot))
  {
    current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
    grad_knot_vals[1,] <- current_grad_val_knot
    xu_vals <- xu

    ## get bounds for knots
    knot_bounds_l <- apply(X = xy, MARGIN = 2, FUN = min)
    knot_bounds_u <- apply(X = xy, MARGIN = 2, FUN = max)
    diffs <- knot_bounds_u - knot_bounds_l
    knot_bounds <- cbind(knot_bounds_l - diffs/10, knot_bounds_u + diffs/10)

    current_xu_trans <- temp_grad_eval$trans_knot
  }
  if(!is.function(dcov_fun_dknot))
  {
    current_grad_val_knot <- 0
  }


  ## run gradient ascent updates until convergence criteria is met or the maxit is reached
  iter <- 1
  # obj_fun_vals[1] <- 1e4 ## this is a line to make the objective fn check work in while loop...

  ## optimization for optim_method = "ga"
  if(optim_method == "ga")
  {
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta, current_grad_val_knot)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {

      ## update the iteration counter
      iter <- iter + 1

      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta, current_grad_val_knot))
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a GA step in transformed space
        if(is.list(dcov_fun_dtheta))
        {
          current_cov_par_trans <- current_cov_par_trans + optim_par$learn_rate * current_grad_val_theta

          ## recover the untransformed parameter values
          for(j in 1:length(cov_par))
          {
            ## transform the transformed parameters back to the original parameter space
            if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
            {
              current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
            }
            else{
              current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                  cov_par = cov_par_start,
                                                                                                                  transform = TRUE)$trans_fun(current_cov_par_trans[j])
            }

          }
          cov_par <- as.list(current_cov_par)

        }

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
          current_xu_trans <- current_xu_trans + optim_par$learn_rate * matrix(data = current_grad_val_knot,
                                                                                   nrow = nrow(xu),
                                                                                   ncol = ncol(xu),
                                                                                   byrow = TRUE)
          xu <- apply(X = current_xu_trans, MARGIN = 1,
                      FUN = dcov_fun_dknot(x1 = 0, x2 = 0, cov_par = cov_par, transform = TRUE, bounds = knot_bounds)$trans_fun,
                      bounds = knot_bounds)

          xu <- matrix(xu, ncol = ncol(xy), byrow = TRUE)

          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }
      }
      if(transform == FALSE)
      {
        current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
        cov_par <- as.list(current_cov_par)

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
          xu <- xu + optim_par$learn_rate* matrix(data = current_grad_val_knot,
                                                nrow = nrow(xu),
                                                ncol = ncol(xu),
                                                byrow = TRUE)
          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }
      }

      ## current objective function value
      Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = as.list(current_cov_par), cov_fun = cov_fun, delta = delta)
      Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(),
                               cov_par = as.list(current_cov_par),
                               cov_fun = cov_fun, delta = delta) - as.list(cov_par)$tau^2 * diag(nrow(xu))

      ## create Z
      Z2 <- solve(a = Sigma22, b = t(Sigma12))
      Z3 <- Sigma12 * t(Z2)
      Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
      Z <- cov_par$sigma^2 + cov_par$tau^2 + delta - Z4

      current_obj_fun <- obj_fun(mu = mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, ...) ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)

      temp_grad_eval <- dlogp_dcov_par(cov_par = as.list(current_cov_par), cov_fun = cov_fun, ## store the results of gradient function
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = knot_opt,
                                       xu = xu, xy = xy,
                                       y = y, ff = NA, mu = mu, delta = delta, transform = transform)
      if(is.list(dcov_fun_dtheta))
      {
        current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
        current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
        grad_vals <- rbind(grad_vals, current_grad_val_theta)
      }

      ## potentially record the gradient of the knot locations
      if(is.function(dcov_fun_dknot))
      {
        current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
        grad_knot_vals <- rbind(grad_knot_vals, current_grad_val_knot)
      }

    }
  }

  ## optimization for optim_method = "adadelta"
  sg2_theta <- 0
  sg2_knot <- 0
  sd2_theta <- 0
  sd2_knot <- 0
  if(optim_method == "adadelta")
  {
    if(is.list(dcov_fun_dtheta))
    {
      ## store the decaying sum of squared gradients
      sg2_theta <- rep(0, times = length(current_cov_par))

      ## store the decaying sum of updates
      sd2_theta <- rep(0, times = length(current_cov_par))

      sign_change_theta <-  rep(0, times = length(current_cov_par))
    }
    if(is.function(dcov_fun_dknot))
    {
      ## store the decaying sum of squared gradients
      sg2_knot <- rep(0, times = nrow(xu)*ncol(xu))

      ## store the decaying sum of updates
      sd2_knot <- rep(0, times = nrow(xu)*ncol(xu))

      sign_change_knot <- rep(0, times = nrow(xu)*ncol(xu))
    }

    ## do the optimization
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta, current_grad_val_knot)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {
      ## update the iteration counter
      iter <- iter + 1
      #
      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta, current_grad_val_knot))
        # print(current_cov_par)
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a adadelta step in transformed space
        if(is.list(dcov_fun_dtheta))
        {
          ## accumulate the gradient
          sg2_theta <- optim_par$decay * sg2_theta + (1 - optim_par$decay) * current_grad_val_theta^2

          ## calculate adapted step size
          # delta_theta <- (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta
          delta_theta <- ((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta

          # delta_theta <- optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))) * current_grad_val_theta

          # if(verbose == TRUE)
          # {
          #   print("adaptive step size theta")
          #   # print(optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))))
          #   print(((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))))
          #   # print(((1/optim_par$eta)^(sign_change_theta)))
          # }

          ## accumulate the step size
          sd2_theta <- optim_par$decay * sd2_theta + (1 - optim_par$decay) * delta_theta^2

          ## update the parameters
          current_cov_par_trans <- current_cov_par_trans + delta_theta
        }

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )

          ## accumulate the gradient
          sg2_knot <- optim_par$decay * sg2_knot + (1 - optim_par$decay) * current_grad_val_knot^2

          ## calculate adapted step size
          # delta_knot <- (sqrt(sd2_knot + rep(optim_par$epsilon, times = length(sd2_knot))) / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot

          delta_knot <- ((1/optim_par$eta)^(sign_change_knot)) * (sqrt(sd2_knot + rep(optim_par$epsilon, times = length(sd2_knot))) / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot
          # delta_knot <- (optim_par$learn_rate / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot

          ## accumulate the step size
          sd2_knot <- optim_par$decay * sd2_knot + (1 - optim_par$decay) * delta_knot^2


          current_xu_trans <- current_xu_trans + matrix(data = delta_knot,
                                                        nrow = nrow(xu),
                                                        ncol = ncol(xu),
                                                        byrow = TRUE)
          xu <- apply(X = current_xu_trans, MARGIN = 1,
                      FUN = dcov_fun_dknot(x1 = 0, x2 = 0, cov_par = cov_par, transform = TRUE, bounds = knot_bounds)$trans_fun,
                      bounds = knot_bounds)

          xu <- matrix(xu, ncol = ncol(xy), byrow = TRUE)

          xu_vals <- abind::abind(xu_vals, xu, along = 3)

        }

        if(is.list(dcov_fun_dtheta))
        {
          ## recover the untransformed parameter values
          for(j in 1:length(cov_par))
          {
            ## transform the transformed parameters back to the original parameter space
            if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
            {
              current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
            }
            else{
              current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                  cov_par = cov_par_start,
                                                                                                                  transform = TRUE)$trans_fun(current_cov_par_trans[j])
            }

          }
          cov_par <- as.list(current_cov_par)

        }
      }
      ## This commented code does gradient ascent... not modified adadelta
      ##  I should probably cut out the transform argument and mandate a transformation to the whole real line
      # if(transform == FALSE)
      # {
      #   current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
      #
      #   if(is.function(dcov_fun_dknot))
      #   {
      #     ## Note that the gradient vector for the knots looks like this
      #     ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
      #     xu <- xu + optim_par$learn_rate * matrix(data = current_grad_val_knot,
      #                                         nrow = nrow(xu),
      #                                         ncol = ncol(xu),
      #                                         byrow = TRUE)
      #     xu_vals <- abind::abind(xu_vals, xu, along = 3)
      #   }
      # }

      ## current objective function value
      if(cov_fun == "ard")
      {
        lnames <- paste("l", 1:ncol(xy), sep = "")
        Sigma12 <- make_cov_mat_ardC(x = xy, x_pred = xu, cov_par = cov_par,
                                     cov_fun = cov_fun,
                                     delta = delta, lnames = lnames)
        Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(),
                                     cov_par = cov_par,
                                     cov_fun = cov_fun,
                                     delta = delta,
                                     lnames = lnames) -
          as.list(cov_par)$tau^2 * diag(nrow(xu))
      }
      else{
        lnames <- character()
        Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
        Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(),
                                 cov_par = cov_par,
                                 cov_fun = cov_fun, delta = delta) - as.list(cov_par)$tau^2 * diag(nrow(xu))
      }

      ## create Z
      Z2 <- solve(a = Sigma22, b = t(Sigma12))
      Z3 <- Sigma12 * t(Z2)
      Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
      Z <- cov_par$sigma^2 + cov_par$tau^2 + delta - Z4

      current_obj_fun <- obj_fun(mu = mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, ...) ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)
      temp_grad_eval <- dlogp_dcov_par(cov_par = as.list(current_cov_par), cov_fun = cov_fun, ## store the results of gradient function
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = knot_opt,
                                       xu = xu, xy = xy,
                                       y = y, ff = NA, mu = mu, delta = delta, transform = transform)


      if(is.list(dcov_fun_dtheta))
      {
        ## get the sign change and add it to the last one
        sign_change_theta <- (optim_par$decay * sign_change_theta +
                                (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))/2)
        # sign_change_theta <- sign_change_theta +
        #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))

        current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
        current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
        grad_vals <- rbind(grad_vals, current_grad_val_theta)
      }

      ## potentially record the gradient of the knot locations
      if(is.function(dcov_fun_dknot))
      {
        ## get the sign change and add it to the last one
        sign_change_knot <- (optim_par$decay * sign_change_knot +
                               (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$knot_gradient)) - sign(current_grad_val_knot)))

        # sign_change_knot <- sign_change_knot +
        #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$knot_gradient)) - sign(current_grad_val_knot))

        current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
        grad_knot_vals <- rbind(grad_knot_vals, current_grad_val_knot)
      }

    }
  }

  ## get the posterior distribution of u
  ## Compute Z^(-1) Sigma12
  ZSig12 <- (1/Z) * Sigma12

  ## compute  mu + Sigma21 %*% Sigma11_inv %*% (fhat - mu)
  ##    the posterior mean of u, the GP at the knot location
  # R1 <- chol(Sigma22)
  R1 <- chol(Sigma22 + t(Sigma12) %*% ZSig12)
  # R2 <- chol(solve(Sigma22) + t(Sigma12) %*% ZSig12)

  # u_mean <- muu + as.numeric(t(ZSig12) %*% (ff - mu)) -
  #   as.numeric( t(Sigma12) %*% (ZSig12 %*%
  #   solve(a = Sigma22 + t(Sigma12) %*% ZSig12, b = t(ZSig12) %*% (ff - mu))) )
  u_mean <- muu + as.numeric(t(ZSig12) %*% (y - mu)) -
    as.numeric( t(Sigma12) %*% (ZSig12 %*%
                                  solve( a = R1, b = solve(a = t(R1), b = t(ZSig12) %*% (y - mu)))) )

  ## compute the posterior variance of u, the GP at the knot location
  u_var <- Sigma22 - t(Sigma12) %*% ZSig12 +
    (t(ZSig12) %*% t(solve(a = t(R1), b = t(Sigma12))) %*% solve(a = t(R1), b = t(Sigma12)) %*% ZSig12)

  if(is.function(dcov_fun_dknot))
  {
    return(list("cov_par" = as.list(current_cov_par),
                "cov_fun" = cov_fun,
                "xu" = xu,
                "xy" = xy,
                "mu" = mu,
                "muu" = muu,
                "u_mean" = u_mean,
                "u_var" = u_var,
                "iter" = iter,
                "obj_fun" = obj_fun_vals,
                "grad" = grad_vals,
                "knot_grad" = grad_knot_vals,
                "knot_history" = xu_vals,
                "cov_par_history" = cov_par_vals))
  }
  if(!is.function(dcov_fun_dknot))
  {
    return(list("cov_par" = as.list(current_cov_par),
                "cov_fun" = cov_fun,
                "xu" = xu,
                "xy" = xy,
                "mu" = mu,
                "muu" = muu,
                "u_mean" = u_mean,
                "u_var" = u_var,
                "cov_fun" = cov_fun,
                "iter" = iter,
                "obj_fun" = obj_fun_vals,
                "grad" = grad_vals,
                "knot_grad" = 0,
                "knot_history" = xu,
                "cov_par_history" = cov_par_vals))
  }

}


## Gradient ascent algorithm for optimizing the full GP with Gaussian data wrt the covariance parameters
## depends on Matrix and glmnet for sparse matrix operations
## depends on dlogpfull_dcov_par to compute gradients
#'@export
norm_grad_ascent_full <- function(cov_par_start,
                             cov_fun,
                             dcov_fun_dtheta,
                             xy,
                             y,
                             mu = NA,
                             transform = TRUE,
                             obj_fun,
                             opt = list(),
                             verbose = FALSE,
                             ...)
{
  ## cov_par_start is a list of the starting values of the covariance function parameters
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is a string specifying the covariance function ("sqexp" or "exp")
  ## xy are the observed data locations (matrix where rows are observations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
  ## transform is a logical argument that says whether to optimize on scales such that paramters are unconstrained
  ## obj_fun is the objective function
  ## learn_rate is the learning rate (step size) for the algorithm
  ## maxit is the maximum number of iterations before cutting it off
  ## tol is the absolute tolerance level which dictates
  ##      how small a difference the algorithm is allowed to make before stopping
  ##      to the objective function as well as how close the gradient is allowed to be to zero
  ## ... are argument to be passed to the functions taking derivatives of log(p(y|lambda)) wrt ff

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = Inf, "delta" = 1e-6)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: you supplied an argument to opt() not utilized by this function.")
        next
      }
      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  optim_method = opt_master$optim_method
  optim_par = list("decay" = opt_master$decay, "epsilon" = opt_master$epsilon, "eta" = opt_master$eta, "learn_rate" = opt_master$learn_rate)
  maxit = opt_master$maxit
  obj_tol = opt_master$obj_tol
  grad_tol = opt_master$grad_tol
  delta <- opt_master$delta

  if(!is.numeric(mu))
  {
    mu <- rep(mean(y), times = length(y))
  }

  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
  }
  else{
    lnames <- character()
  }

  ## make sure y is numeric
  y <- as.numeric(y)

  ## store objective function evaluations, gradient values, and parameter values
  ## initialize vector of objective function values
  current_cov_par <- unlist(cov_par_start) ## vector of covariance parameter values (without inducing points)
  cov_par_vals <- matrix(ncol = length(current_cov_par)) ## store the covariance parameter values
  cov_par <- as.list(current_cov_par)

  obj_fun_vals <-  numeric() ## store the objective function values
  current_obj_fun <- NA ## current objective function value

  grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
  current_grad_val_theta <- NA ## current gradient value

  ## create first covariance matrix
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
    Sigma11 <- make_cov_mat_ardC(x = xy, x_pred = matrix(),
                                 cov_par = cov_par,
                                 cov_fun = cov_fun,
                                 delta = delta,
                                 lnames = lnames)

  }
  else{
    lnames <- character()
    Sigma11 <- make_cov_matC(x = xy, x_pred = matrix(), cov_par = cov_par, cov_fun = cov_fun, delta = delta)
  }

  ## current objective function value
  current_obj_fun <- obj_fun(mu = mu, Sigma = Sigma11, y = y, ...) ## store the first objective function value
  obj_fun_vals[1] <- current_obj_fun
  cov_par_vals[1,] <- current_cov_par
  temp_grad_eval <- dlogp_dcov_par_full(cov_par = cov_par, cov_fun = cov_fun, ## store the results of gradient function
                                   dcov_fun_dtheta = dcov_fun_dtheta,
                                   xy = xy,
                                   y = y, mu = mu, delta = delta, transform = transform)

  current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
  grad_vals[1,] <- current_grad_val_theta
  current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters

  ## run gradient ascent updates until convergence criteria is met or the maxit is reached
  iter <- 1
  # obj_fun_vals[1] <- 1e4 ## this is a line to make the objective fn check work in while loop...

  ## optimization for optim_method = "ga"
  if(optim_method == "ga")
  {
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {

      ## update the iteration counter
      iter <- iter + 1

      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta))
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a GA step in transformed space
        current_cov_par_trans <- current_cov_par_trans + optim_par$learn_rate * current_grad_val_theta

        ## recover the untransformed parameter values
        for(j in 1:length(cov_par))
        {
          ## transform the transformed parameters back to the original parameter space
          if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
          {
            current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
          }
          else{
            current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                cov_par = cov_par_start,
                                                                                                                transform = TRUE)$trans_fun(current_cov_par_trans[j])
          }

        }
        cov_par <- as.list(current_cov_par)

      }
      if(transform == FALSE)
      {
        current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
        cov_par <- as.list(current_cov_par)

      }

      ## current objective function value
      Sigma11 <- make_cov_matC(x = xy, x_pred = matrix(), cov_par = cov_par, cov_fun = cov_fun, delta = delta)

      current_obj_fun <- obj_fun(mu = mu, Sigma = Sigma11, y = y, ...) ## store the objective function value
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)

      temp_grad_eval <- dlogp_dcov_par_full(cov_par = cov_par, cov_fun = cov_fun, ## store the results of gradient function
                                            dcov_fun_dtheta = dcov_fun_dtheta,
                                            xy = xy,
                                            y = y, mu = mu, delta = delta, transform = transform)
      current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
      current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
      grad_vals <- rbind(grad_vals, current_grad_val_theta)

    }
  }

  ## optimization for optim_method = "adadelta"
  sg2_theta <- 0
  sg2_knot <- 0
  sd2_theta <- 0
  sd2_knot <- 0
  if(optim_method == "adadelta")
  {
    if(is.list(dcov_fun_dtheta))
    {
      ## store the decaying sum of squared gradients
      sg2_theta <- rep(0, times = length(current_cov_par))

      ## store the decaying sum of updates
      sd2_theta <- rep(0, times = length(current_cov_par))

      sign_change_theta <-  rep(0, times = length(current_cov_par))
    }

    ## do the optimization
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {
      ## update the iteration counter
      iter <- iter + 1
      #
      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta))
        # print(current_cov_par)
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a adadelta step in transformed space
        ## accumulate the gradient
        sg2_theta <- optim_par$decay * sg2_theta + (1 - optim_par$decay) * current_grad_val_theta^2

        ## calculate adapted step size
        # delta_theta <- (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta
        delta_theta <- ((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta

        # delta_theta <- optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))) * current_grad_val_theta

        # if(verbose == TRUE)
        # {
        #   print("adaptive step size theta")
        #   # print(optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))))
        #   print(((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))))
        #   # print(((1/optim_par$eta)^(sign_change_theta)))
        # }

        ## accumulate the step size
        sd2_theta <- optim_par$decay * sd2_theta + (1 - optim_par$decay) * delta_theta^2

        ## update the parameters
        current_cov_par_trans <- current_cov_par_trans + delta_theta

        ## recover the untransformed parameter values
        for(j in 1:length(cov_par))
        {
          ## transform the transformed parameters back to the original parameter space
          if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
          {
            current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
          }
          else{
            current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                cov_par = cov_par_start,
                                                                                                                transform = TRUE)$trans_fun(current_cov_par_trans[j])
          }

        }
        cov_par <- as.list(current_cov_par)

      }

      ## current objective function value
      if(cov_fun == "ard")
      {
        lnames <- paste("l", 1:ncol(xy), sep = "")
        Sigma11 <- make_cov_mat_ardC(x = xy, x_pred = matrix(),
                                     cov_par = cov_par,
                                     cov_fun = cov_fun,
                                     delta = delta, lnames = lnames)
      }
      else{
        lnames <- character()
        Sigma11 <- make_cov_matC(x = xy, x_pred = matrix(), cov_par = cov_par, cov_fun = cov_fun, delta = delta)
      }

      current_obj_fun <- obj_fun(mu = mu, Sigma = Sigma11, y = y, ...) ## store the objective function value
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)

      temp_grad_eval <- dlogp_dcov_par_full(cov_par = cov_par, cov_fun = cov_fun, ## store the results of gradient function
                                            dcov_fun_dtheta = dcov_fun_dtheta,
                                            xy = xy,
                                            y = y, mu = mu, delta = delta, transform = transform)

      ## get the sign change and add it to the last one
      sign_change_theta <- (optim_par$decay * sign_change_theta +
                              (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))/2)
      # sign_change_theta <- sign_change_theta +
      #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))

      current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
      current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
      grad_vals <- rbind(grad_vals, current_grad_val_theta)

    }
  }

  return(list("cov_par" = as.list(current_cov_par),
              "cov_fun" = cov_fun,
              "xy" = xy,
              "y" = y,
              "mu" = mu,
              "iter" = iter,
              "obj_fun" = obj_fun_vals,
              "grad" = grad_vals,
              "cov_par_history" = cov_par_vals))

}

## SOR gradient ascent
## Gradient ascent algorithm for optimizing the likelihood wrt the covariance parameters (and eventually inducing points)
## depends on Matrix and glmnet for sparse matrix operations and abind::abind for keeping track of xu values at each iteration
## depends on dlogq_dcov_par to compute gradients of Laplace Approximation
#'@export
norm_grad_ascent_sor <- function(cov_par_start,
                             cov_fun,
                             dcov_fun_dtheta,
                             dcov_fun_dknot,
                             knot_opt,
                             xu,
                             xy,
                             y,
                             mu = NA,
                             muu = NA,
                             transform = TRUE,
                             obj_fun,
                             opt = list(),
                             verbose = FALSE,
                             ...)
{
  ## cov_par_start is a list of the starting values of the covariance function parameters
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is a string specifying the covariance function ("sqexp" or "exp")
  ## xy are the observed data locations (matrix where rows are observations)
  ## knot_opt is a vector giving the knots over which to optimize, other knots are held constant
  ## xu are the unobserved knot locations (matrix where rows are knot locations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## muu is the mean of the GP at each of the knot locations
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

  ## jiggle xu to make sure that knot locations aren't exactly at data locations
  ## only necessary so that the covariance function doesn't add a nugget
  ## between two nuggetless function values
  for(i in 1:nrow(xu))
  {
    if(any(
      apply(X = xy, MARGIN = 1, FUN = function(x,y){isTRUE(all.equal(target = y, current = x))}, y = xu[i,]) == TRUE
    )
    )
    {
      xu[i,] <- xu[i,] + rnorm(n = ncol(xu),sd = 1e-6)
    }
  }

  ## extract names of optional arguments from opt
  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = Inf, "delta" = 1e-6)

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        # print("Warning: you supplied an argument to opt() not utilized by this function.")
        next
      }
      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  optim_method = opt_master$optim_method
  optim_par = list("decay" = opt_master$decay, "epsilon" = opt_master$epsilon, "eta" = opt_master$eta, "learn_rate" = opt_master$learn_rate)
  maxit = opt_master$maxit
  obj_tol = opt_master$obj_tol
  grad_tol = opt_master$grad_tol
  delta <- opt_master$delta

  if(!is.numeric(mu))
  {
    mu <- rep(mean(y), times = length(y))
  }
  if(!is.numeric(muu))
  {
    muu <- rep(mean(y), times = nrow(xu))
  }

  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
  }
  else{
    lnames <- character()
  }

  ## make sure y is numeric
  y <- as.numeric(y)

  ## store objective function evaluations, gradient values, and parameter values
  ## initialize vector of objective function values
  current_cov_par <- unlist(cov_par_start) ## vector of covariance parameter values (without inducing points)
  cov_par_vals <- matrix(ncol = length(current_cov_par)) ## store the covariance parameter values
  cov_par <- as.list(current_cov_par)

  obj_fun_vals <-  numeric() ## store the objective function values
  current_obj_fun <- NA ## current objective function value

  if(is.list(dcov_fun_dtheta))
  {
    grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
    current_grad_val_theta <- NA ## current gradient value
  }
  if(!is.list(dcov_fun_dtheta))
  {
    current_grad_val_theta <- 0 ## current gradient value
    grad_vals <- matrix(ncol = length(current_cov_par)) ## store gradient values
  }

  grad_knot_vals <- matrix(ncol = nrow(xu) * ncol(xu)) ## store gradient values
  current_grad_val_knot <- NA ## current gradient value

  ## current objective function value
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
    Sigma12 <- make_cov_mat_ardC(x = xy, x_pred = xu, cov_par = cov_par,
                                 cov_fun = cov_fun,
                                 delta = delta, lnames = lnames)
    Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(),
                                 cov_par = cov_par,
                                 cov_fun = cov_fun,
                                 delta = delta,
                                 lnames = lnames) -
      as.list(cov_par)$tau^2 * diag(nrow(xu))
  }
  else{
    lnames <- character()
    Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
    Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(),
                             cov_par = cov_par,
                             cov_fun = cov_fun, delta = delta) - as.list(cov_par)$tau^2 * diag(nrow(xu))
  }

  ## create Z
  # Z2 <- solve(a = Sigma22, b = t(Sigma12))
  Z <- rep(cov_par$tau^2 + delta , times = nrow(Sigma12))

  current_obj_fun <- obj_fun(mu = mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, ...) ## store the first objective function value log(q(y|theta, xu, fmax))
  obj_fun_vals[1] <- current_obj_fun
  cov_par_vals[1,] <- current_cov_par
  temp_grad_eval <- dlogp_dcov_par_sor(cov_par = cov_par, cov_fun = cov_fun, ## store the results of gradient function
                                   dcov_fun_dtheta = dcov_fun_dtheta,
                                   dcov_fun_dknot = dcov_fun_dknot,
                                   knot_opt = knot_opt,
                                   xu = xu, xy = xy,
                                   y = y, ff = NA, mu = mu, delta = delta, transform = transform)

  if(is.list(dcov_fun_dtheta))
  {
    current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
    grad_vals[1,] <- current_grad_val_theta
  }
  current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters


  ## potentially record the gradient of the knot locations
  if(is.function(dcov_fun_dknot))
  {
    current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
    grad_knot_vals[1,] <- current_grad_val_knot
    xu_vals <- xu

    ## get bounds for knots
    knot_bounds_l <- apply(X = xy, MARGIN = 2, FUN = min)
    knot_bounds_u <- apply(X = xy, MARGIN = 2, FUN = max)
    diffs <- knot_bounds_u - knot_bounds_l
    knot_bounds <- cbind(knot_bounds_l - diffs/10, knot_bounds_u + diffs/10)

    current_xu_trans <- temp_grad_eval$trans_knot
  }
  if(!is.function(dcov_fun_dknot))
  {
    current_grad_val_knot <- 0
  }


  ## run gradient ascent updates until convergence criteria is met or the maxit is reached
  iter <- 1
  # obj_fun_vals[1] <- 1e4 ## this is a line to make the objective fn check work in while loop...

  ## optimization for optim_method = "ga"
  if(optim_method == "ga")
  {
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta, current_grad_val_knot)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {

      ## update the iteration counter
      iter <- iter + 1

      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta, current_grad_val_knot))
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a GA step in transformed space
        if(is.list(dcov_fun_dtheta))
        {
          current_cov_par_trans <- current_cov_par_trans + optim_par$learn_rate * current_grad_val_theta

          ## recover the untransformed parameter values
          for(j in 1:length(cov_par))
          {
            ## transform the transformed parameters back to the original parameter space
            if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
            {
              current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
            }
            else{
              current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                  cov_par = cov_par_start,
                                                                                                                  transform = TRUE)$trans_fun(current_cov_par_trans[j])
            }

          }
          cov_par <- as.list(current_cov_par)

        }

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
          current_xu_trans <- current_xu_trans + optim_par$learn_rate * matrix(data = current_grad_val_knot,
                                                                               nrow = nrow(xu),
                                                                               ncol = ncol(xu),
                                                                               byrow = TRUE)
          xu <- apply(X = current_xu_trans, MARGIN = 1,
                      FUN = dcov_fun_dknot(x1 = 0, x2 = 0, cov_par = cov_par, transform = TRUE, bounds = knot_bounds)$trans_fun,
                      bounds = knot_bounds)

          xu <- matrix(xu, ncol = ncol(xy), byrow = TRUE)

          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }
      }
      if(transform == FALSE)
      {
        current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
        cov_par <- as.list(current_cov_par)

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
          xu <- xu + optim_par$learn_rate* matrix(data = current_grad_val_knot,
                                                  nrow = nrow(xu),
                                                  ncol = ncol(xu),
                                                  byrow = TRUE)
          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }
      }

      ## current objective function value
      Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = as.list(current_cov_par), cov_fun = cov_fun, delta = delta)
      Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(),
                               cov_par = as.list(current_cov_par),
                               cov_fun = cov_fun, delta = delta) - as.list(cov_par)$tau^2 * diag(nrow(xu))

      ## create Z
      Z2 <- solve(a = Sigma22, b = t(Sigma12))
      Z <- numeric()
      Z <- rep(cov_par$tau^2 + delta, times = nrow(Sigma12))

      current_obj_fun <- obj_fun(mu = mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, ...) ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)

      temp_grad_eval <- dlogp_dcov_par_sor(cov_par = as.list(current_cov_par), cov_fun = cov_fun, ## store the results of gradient function
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = knot_opt,
                                       xu = xu, xy = xy,
                                       y = y, ff = NA, mu = mu, delta = delta, transform = transform)
      if(is.list(dcov_fun_dtheta))
      {
        current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
        current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
        grad_vals <- rbind(grad_vals, current_grad_val_theta)
      }

      ## potentially record the gradient of the knot locations
      if(is.function(dcov_fun_dknot))
      {
        current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
        grad_knot_vals <- rbind(grad_knot_vals, current_grad_val_knot)
      }

    }
  }

  ## optimization for optim_method = "adadelta"
  sg2_theta <- 0
  sg2_knot <- 0
  sd2_theta <- 0
  sd2_knot <- 0
  if(optim_method == "adadelta")
  {
    if(is.list(dcov_fun_dtheta))
    {
      ## store the decaying sum of squared gradients
      sg2_theta <- rep(0, times = length(current_cov_par))

      ## store the decaying sum of updates
      sd2_theta <- rep(0, times = length(current_cov_par))

      sign_change_theta <-  rep(0, times = length(current_cov_par))
    }
    if(is.function(dcov_fun_dknot))
    {
      ## store the decaying sum of squared gradients
      sg2_knot <- rep(0, times = nrow(xu)*ncol(xu))

      ## store the decaying sum of updates
      sd2_knot <- rep(0, times = nrow(xu)*ncol(xu))

      sign_change_knot <- rep(0, times = nrow(xu)*ncol(xu))
    }

    ## do the optimization
    while(iter < maxit &&
          (any(abs(c(current_grad_val_theta, current_grad_val_knot)) > grad_tol) ||
           ifelse(iter > 1, yes = abs(current_obj_fun - obj_fun_vals[iter - 1]) > obj_tol, no = TRUE)))
    {
      ## update the iteration counter
      iter <- iter + 1
      #
      if(verbose == TRUE)
      {
        print(paste("iteration ", iter, sep = ""))
        print(c(current_grad_val_theta, current_grad_val_knot))
        # print(current_cov_par)
        # print(current_obj_fun)
      }
      if(transform == TRUE)
      {
        ## take a adadelta step in transformed space
        if(is.list(dcov_fun_dtheta))
        {
          ## accumulate the gradient
          sg2_theta <- optim_par$decay * sg2_theta + (1 - optim_par$decay) * current_grad_val_theta^2

          ## calculate adapted step size
          # delta_theta <- (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta
          delta_theta <- ((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))) * current_grad_val_theta

          # delta_theta <- optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))) * current_grad_val_theta

          # if(verbose == TRUE)
          # {
          #   print("adaptive step size theta")
          #   # print(optim_par$learn_rate / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta))))
          #   print(((1/optim_par$eta)^(sign_change_theta)) * (sqrt(sd2_theta + rep(optim_par$epsilon, times = length(sd2_theta)) ) / sqrt(sg2_theta + rep(optim_par$epsilon, times = length(sd2_theta)))))
          #   # print(((1/optim_par$eta)^(sign_change_theta)))
          # }

          ## accumulate the step size
          sd2_theta <- optim_par$decay * sd2_theta + (1 - optim_par$decay) * delta_theta^2

          ## update the parameters
          current_cov_par_trans <- current_cov_par_trans + delta_theta
        }

        if(is.function(dcov_fun_dknot))
        {
          ## Note that the gradient vector for the knots looks like this
          ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )

          ## accumulate the gradient
          sg2_knot <- optim_par$decay * sg2_knot + (1 - optim_par$decay) * current_grad_val_knot^2

          ## calculate adapted step size
          # delta_knot <- (sqrt(sd2_knot + rep(optim_par$epsilon, times = length(sd2_knot))) / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot

          delta_knot <- ((1/optim_par$eta)^(sign_change_knot)) * (sqrt(sd2_knot + rep(optim_par$epsilon, times = length(sd2_knot))) / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot
          # delta_knot <- (optim_par$learn_rate / sqrt(sg2_knot + rep(optim_par$epsilon, times = length(sd2_knot)))) * current_grad_val_knot

          ## accumulate the step size
          sd2_knot <- optim_par$decay * sd2_knot + (1 - optim_par$decay) * delta_knot^2


          current_xu_trans <- current_xu_trans + matrix(data = delta_knot,
                                                        nrow = nrow(xu),
                                                        ncol = ncol(xu),
                                                        byrow = TRUE)
          xu <- apply(X = current_xu_trans, MARGIN = 1,
                      FUN = dcov_fun_dknot(x1 = 0, x2 = 0, cov_par = cov_par, transform = TRUE, bounds = knot_bounds)$trans_fun,
                      bounds = knot_bounds)

          xu <- matrix(xu, ncol = ncol(xy), byrow = TRUE)

          xu_vals <- abind::abind(xu_vals, xu, along = 3)

        }

        if(is.list(dcov_fun_dtheta))
        {
          ## recover the untransformed parameter values
          for(j in 1:length(cov_par))
          {
            ## transform the transformed parameters back to the original parameter space
            if(cov_fun == "ard" & names(cov_par)[j] %in% lnames)
            {
              current_cov_par[j] <- real_to_pos(x = current_cov_par_trans[j])
            }
            else{
              current_cov_par[j] <- dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == names(current_cov_par)[j])]](x1 = 0, x2 = 2,
                                                                                                                  cov_par = cov_par_start,
                                                                                                                  transform = TRUE)$trans_fun(current_cov_par_trans[j])
            }

          }
          cov_par <- as.list(current_cov_par)

        }
      }
      ## This commented code does gradient ascent... not modified adadelta
      ##  I should probably cut out the transform argument and mandate a transformation to the whole real line
      # if(transform == FALSE)
      # {
      #   current_cov_par <- current_cov_par + optim_par$learn_rate * current_grad_val_theta
      #
      #   if(is.function(dcov_fun_dknot))
      #   {
      #     ## Note that the gradient vector for the knots looks like this
      #     ## c( xu[1,1], xu[1,2], ..., xu[1,d], xu[2,1], ..., xu[m,d] )
      #     xu <- xu + optim_par$learn_rate * matrix(data = current_grad_val_knot,
      #                                         nrow = nrow(xu),
      #                                         ncol = ncol(xu),
      #                                         byrow = TRUE)
      #     xu_vals <- abind::abind(xu_vals, xu, along = 3)
      #   }
      # }

      ## current objective function value
      if(cov_fun == "ard")
      {
        lnames <- paste("l", 1:ncol(xy), sep = "")
        Sigma12 <- make_cov_mat_ardC(x = xy, x_pred = xu, cov_par = cov_par,
                                     cov_fun = cov_fun,
                                     delta = delta, lnames = lnames)
        Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(),
                                     cov_par = cov_par,
                                     cov_fun = cov_fun,
                                     delta = delta,
                                     lnames = lnames) -
          as.list(cov_par)$tau^2 * diag(nrow(xu))
      }
      else{
        lnames <- character()
        Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
        Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(),
                                 cov_par = cov_par,
                                 cov_fun = cov_fun, delta = delta) - as.list(cov_par)$tau^2 * diag(nrow(xu))
      }

      ## create Z
      Z2 <- solve(a = Sigma22, b = t(Sigma12))
      Z <- numeric()
      Z <- rep(cov_par$tau^2 + delta, times = nrow(Sigma12))

      current_obj_fun <- obj_fun(mu = mu, Z = Z, Sigma12 = Sigma12, Sigma22 = Sigma22, y = y, ...) ## store the first objective function value log(q(y|theta, xu, fmax))
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)
      temp_grad_eval <- dlogp_dcov_par_sor(cov_par = as.list(current_cov_par), cov_fun = cov_fun, ## store the results of gradient function
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = knot_opt,
                                       xu = xu, xy = xy,
                                       y = y, ff = NA, mu = mu, delta = delta, transform = transform)


      if(is.list(dcov_fun_dtheta))
      {
        ## get the sign change and add it to the last one
        sign_change_theta <- (optim_par$decay * sign_change_theta +
                                (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))/2)
        # sign_change_theta <- sign_change_theta +
        #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$gradient)) - sign(current_grad_val_theta))

        current_grad_val_theta <- c(temp_grad_eval$gradient) ## store the current gradient
        current_cov_par_trans <- unlist(temp_grad_eval$trans_par) ## store the current value of the transformed covariance parameters
        grad_vals <- rbind(grad_vals, current_grad_val_theta)
      }

      ## potentially record the gradient of the knot locations
      if(is.function(dcov_fun_dknot))
      {
        ## get the sign change and add it to the last one
        sign_change_knot <- (optim_par$decay * sign_change_knot +
                               (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$knot_gradient)) - sign(current_grad_val_knot)))

        # sign_change_knot <- sign_change_knot +
        #   (1 - optim_par$decay) * abs(sign(c(temp_grad_eval$knot_gradient)) - sign(current_grad_val_knot))

        current_grad_val_knot <- c(temp_grad_eval$knot_gradient) ## store the current gradient
        grad_knot_vals <- rbind(grad_knot_vals, current_grad_val_knot)
      }

    }
  }

  ## get the posterior distribution of u
  ## Compute Z^(-1) Sigma12
  ZSig12 <- (1/Z) * Sigma12

  ## compute  mu + Sigma21 %*% Sigma11_inv %*% (fhat - mu)
  ##    the posterior mean of u, the GP at the knot location
  # R1 <- chol(Sigma22)
  R1 <- chol(Sigma22 + t(Sigma12) %*% ZSig12)
  # R2 <- chol(solve(Sigma22) + t(Sigma12) %*% ZSig12)

  # u_mean <- muu + as.numeric(t(ZSig12) %*% (ff - mu)) -
  #   as.numeric( t(Sigma12) %*% (ZSig12 %*%
  #   solve(a = Sigma22 + t(Sigma12) %*% ZSig12, b = t(ZSig12) %*% (ff - mu))) )
  u_mean <- muu + as.numeric(t(ZSig12) %*% (y - mu)) -
    as.numeric( t(Sigma12) %*% (ZSig12 %*%
                                  solve( a = R1, b = solve(a = t(R1), b = t(ZSig12) %*% (y - mu)))) )

  ## compute the posterior variance of u, the GP at the knot location
  u_var <- Sigma22 - t(Sigma12) %*% ZSig12 +
    (t(ZSig12) %*% t(solve(a = t(R1), b = t(Sigma12))) %*% solve(a = t(R1), b = t(Sigma12)) %*% ZSig12)

  if(is.function(dcov_fun_dknot))
  {
    return(list("cov_par" = as.list(current_cov_par),
                "cov_fun" = cov_fun,
                "xu" = xu,
                "xy" = xy,
                "mu" = mu,
                "muu" = muu,
                "u_mean" = u_mean,
                "u_var" = u_var,
                "iter" = iter,
                "obj_fun" = obj_fun_vals,
                "grad" = grad_vals,
                "knot_grad" = grad_knot_vals,
                "knot_history" = xu_vals,
                "cov_par_history" = cov_par_vals))
  }
  if(!is.function(dcov_fun_dknot))
  {
    return(list("cov_par" = as.list(current_cov_par),
                "cov_fun" = cov_fun,
                "xu" = xu,
                "xy" = xy,
                "mu" = mu,
                "muu" = muu,
                "u_mean" = u_mean,
                "u_var" = u_var,
                "cov_fun" = cov_fun,
                "iter" = iter,
                "obj_fun" = obj_fun_vals,
                "grad" = grad_vals,
                "knot_grad" = 0,
                "knot_history" = xu,
                "cov_par_history" = cov_par_vals))
  }

}
