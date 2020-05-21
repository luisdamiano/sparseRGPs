#' Train a GP via type-II maximum likelihood.
#'
#' @description Train a full or sparse GP via type-II maximum likelihood.
#' @param y Vector of the observed response values.
#' @param xy Matrix of the observed input/covariate values. Rows correspond to
#' elements of y.
#' @param cov_fun Chacter string specifying the GP covariance function.
#'  Either "sqexp" for squared exponential, or "ard" for automatic
#'  relevance determination are currently supported.
#' @param cov_par_start Named list of initial values for the covariance
#' parameters. Always requires "sigma" and "tau" which are the standard
#' deviations for the latent GP and the noise, respectively. In the case
#' of the squared exponential covariance function, "l" must also be specified.
#' In the case of the ARD covariance function, "l1", "l2", ..., "ld"
#' must be specified where d is the dimension of the input space/number of
#' columns of xy.
#' @param mu Vector of the same length of y specifying the mean of the GP.
#' @param family Character string either "gaussian", "bernoulli", or "poisson"
#' specifying the type of data that y is.
#' @param nugget Logical value indicating whether to estimate the nugget alongside
#' other covariance parameters.
#' @param sparse Logical value indicating whether to fit a sparse or full GP.
#' @param xu_opt Character string denoting how to select knots in the case of
#' a sparse model. "fixed" fixes the knots to initial values. "random"
#' uses the OAT knot selection algorithm with a best of random subset proposal.
#' "oat" uses Bayesian optimization to propose a knot. "simultaneous" simultaneously
#' optimizes knots alongside covariance parameters.
#' @param xu Matrix of initial knots. Number of columns should match that of xy.
#' @param muu Vector of the marginal mean of the GP at the knots.
#' @param vi Logical value indicating whether to fit the sparse model using
#' variational inference. This only works for Gaussian data.
#' @param opt A named list of additional options to control the optimization
#' algorithm and the knot selection procedure. The only arguments that
#' will likely need to be set are "delta", "maxknot", and potentially
#' "TTmax". All others should typically work well with default values.
#'
#'     decay: Controls how quickly previous gradients are forgotten. Values
#'     equal to zero recovers gradient ascent, i.e. previous gradients are
#'     forgotten entirely. Values close to 1 indicate that the current gradient
#'     has little influence over the step direction. Defaults to 0.95.
#'
#'     epsilon: Small, positive value partially controlling the step size.
#'     Defaults to 1e-6.
#'     eta: Additional parameter >= 1 added to Adadelta that shrinks step sizes
#'     when gradients oscillate between negative and positive values.
#'     If equal to 1, recovers Adadelta. Defaults to 1e3.
#'
#'     maxit: Maximum number of gradient ascent iterations.
#'     Defaults to 1000.
#'
#'     obj_tol: Objective function tolerance. Convergence declared when
#'     change in the objective function falls below this value.
#'     Defaults to 1e-3.
#'
#'     grad_tol: Gradient tolerance. Additional tolerance parameter
#'     controlling convergence of gradient ascent. Convergence declared if
#'     objective tolerance is met AND all absolute gradients fall below
#'     a threshold. Defaults to Inf so that only obj_tol controls convergence.
#'
#'     maxit_nr: Maximum number of Newton-Raphson iterations if data is non-Gaussian.
#'
#'     delta: Small, positive quantity added to the marginal GP variances to
#'     ensure well-conditioned matrices. Defaults to 1e-6. Numerical problems
#'     can occur when sigma^2 / (tau^2 + delta) > 1e6.
#'
#'     tol_nr: Objective function and gradient tolerance level for the
#'     Newton-Raphson algorithm. Defaults to 1e-5. Small tolerance levels here
#'     are recommended, but can sometimes slow down training.
#'
#'     maxknot: Maximum number of allowable knots
#'
#'     cov_diff: Logical value indicating whether to treat the covariance
#'     function as if it is not differentiable so that each added knot
#'     is not optimized continuously
#'
#'     chooseK: Logical value indicating whether to select the total number
#'     of knots or continue adding knots until maxknot is reached
#'
#'     TTmax: Maximum number of candidate knot proposals to test in the OAT
#'     algorithm
#'
#'     TTmin: Minimum number of candidate knot proposals to test in the OAT
#'     algorithm
#'
#'     ego_cov_par: Named list of initial covariance parameters values for
#'     the meta GP
#'     in the case that Bayesian optimization is used for knot proposals.
#'     Values should not need to be changed from the default.
#' @param verbose Logical value indicating whether to print the iteration and
#' the current covariance parameter/knot gradient value.
#' @param file_path Character string denoting the path to the file where you
#' would like to save the trained model.
#' @param ... Additional, model dependent arguments. Currently the only use
#' for this is with Poisson response values. In this case, the user should
#' specify the variable m, which is a vector of length equal to y. This
#' is a part of the mean of the Poisson where E(Y) = m*exp(f(x)) and f(x)
#' is the value of the latent GP at x.
#' @return Returns a list of values corresponding to the fitted GP. These include:
#'
#'     sparse: logical value indicating whether the model is sparse or not
#'
#'     family: Character string indicating the conditional distribution of the
#'     data given the latent function
#'
#'     delta: A number corresponding to the value of delta used in the fitted
#'     model to stabilize relevant matrix operations.
#'
#'     xu_init: A matrix of initial knot values in the case that a
#'     sparse model was fit.
#'
#'     results: A list with the following elements:
#'
#'         results$cov_par: A named list with the fitted covariance parameter values
#'
#'         results$cov_fun: Character string indicating the covariance function.
#'
#'         results$xu: The fitted, final knot locations
#'
#'         results$obj_fun: A vector of objective function values for every gradient
#'         ascent step if OAT is not used, or the optimized objective
#'         function after adding each knot if OAT is used.
#'
#'          results$u_mean: The posterior mean of the latent function at the knots.
#'
#'          results$u_var: The posterior variance-covariance matrix of the
#'          latent function values at the knots.
#'
#'          results$muu: The marginal mean of the latent function at the knots.
#'
#'          results$mu: The marginal mean of the latent function at the observed
#'          data locations
#'
#'          results$cov_par_history: Either the covariance parameter values at
#'          every gradient ascent step if OAT is not used, or the optimized
#'          covariance parameter values after adding each knot if OAT
#'          is used.
#'
#'          results$ga_steps: A vector giving the number of gradient ascent steps,
#'          potentially for every added knot if OAT is used.
#'
#'          results$u_post: a list where each element is a list showing the posterior
#'          mean and variance-covariance matrix of the
#'          latent function at the knots after adding each knot
#' @importFrom Rdpack reprompt
#'
#' @export
optimize_gp <- function(y,
                        xy,
                        cov_fun = "sqexp",
                        cov_par_start,
                        mu,
                        family = "gaussian",
                        nugget = TRUE,
                        sparse = FALSE,
                        xu_opt = NULL,
                        xu = NULL,
                        muu = NULL,
                        vi = FALSE,
                        opt = NULL,
                        verbose = FALSE,
                        file_path = NA,
                        ...)
{

  ## extract names of optional arguments from opt
  opt_master = list("optim_method" = "adadelta",
                    "decay" = 0.95, "epsilon" = 1e-6, "learn_rate" = 1e-2, "eta" = 1e3,
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = Inf,
                    "maxit_nr" = 1000, "delta" = 1e-6, "tol_nr" = 1e-5,
                    "maxknot" = 15, "tol_knot" = 1e-3,
                    "ego_cov_par" = list("sigma" = 3, "l" = 1, "tau" = 0),
                    "ego_cov_fun" = "exp",
                    "ego_cov_fun_dtheta" = list("sigma" = dexp_dsigma,
                                                "l" = dexp_dl),
                    "TTmax" = 30, "TTmin" = 10,
                    "cov_diff" = TRUE,
                    "chooseK" = TRUE)

  ## don't let user use basic gradient ascent
  if(opt_master$optim_method == "ga")
  {
    return("Error: Basic gradient ascent no longer supported.")
  }

  ## print error if try to use VI with non gaussian data
  if(vi == TRUE && family != "gaussian")
  {
    return("Error: VI not supported for non-gaussian data.")
  }

  ## potentially change default values in opt_master to user supplied values
  if(length(opt) > 0)
  {
    for(i in 1:length(opt))
    {
      ## if an argument is misspelled or not accepted, throw an error
      if(!any(names(opt_master) == names(opt)[i]))
      {
        print("Warning: you supplied an argument to opt() not utilized by this function.")
      }
      else{
        ind <- which(names(opt_master) == names(opt)[i])
        opt_master[[ind]] <- opt[[i]]
      }

    }
  }

  ## rename opt_master to opt
  opt <- opt_master

  ## check to make sure that 0 < TTmin < TTmax < N - K
  ##    K < N
  if(sparse == TRUE & maxknot >= length(y))
  {
    print("Error: Maximum number of knots must be less than
           the number of data points.")
    return()
  }
  if(!is.null(xu_opt))
  {
    if(xu_opt %in% c("random", "oat") & (
      opt$TTmin <= 0 ||
      opt$TTmax <= opt$TTmin ||
      opt$TTmax >= (length(y) - maxknot)
    )
    ){
      print("Error: Must adhere to 0 < TTmin < TTmax < N - maxknot. Recommend
           setting both TTmax and maxknot < N / 2 or using a full GP.")
      return()
    }
  }


  args <- list(...)
  if(cov_fun == "sqexp")
  {
    if(nugget == TRUE)
    {
      dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl, "tau" = dsqexp_dtau)
    }
    if(nugget == FALSE)
    {
      dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl)
    }
    dcov_fun_dknot <- dsqexp_dx2
  }
  if(cov_fun == "ard")
  {
    if(nugget == TRUE)
    {
      dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma_ard,
                              # "l" = dsqexp_dl,
                              "tau" = dsqexp_dtau)
    }
    if(nugget == FALSE)
    {
      dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma_ard
      )
    }
    dcov_fun_dknot <- dsqexp_dx2_ard
  }

  if(!cov_fun %in% c("sqexp","ard"))
  {
    print("Error: invalid covariance function")
    return()
  }

  if(!is.matrix(xy))
  {
    print("Warning: xy must be a matrix. I'll try to make the conversion.")
    xy <- matrix(data = xy, ncol = 1)
  }

  ## sparse case
  if(sparse == TRUE)
  {
    if(!is.matrix(xu))
    {
      print("Warning: xu must be a matrix. I'll try to make the conversion.")
      xu <- matrix(data = xu, ncol = 1)
    }

    if(length(muu) != nrow(xu))
    {
      print("Warning: length(muu) is not equal to nrow(xu). Setting muu = 0")
      muu <- rep(0, times = nrow(xu))
    }

    ######################################
    ## Gaussian case
    ######################################
    if(family == "gaussian")
    {

      if(vi == TRUE)
      {
        if(xu_opt == "fixed")
        {
          dcov_fun_dknot <- NA
          results <- norm_grad_ascent_vi(cov_par_start = cov_par_start,
                                         cov_fun = cov_fun,
                                         dcov_fun_dtheta = dcov_fun_dtheta,
                                         dcov_fun_dknot = dcov_fun_dknot,
                                         knot_opt = NA, xu = xu,
                                         xy = xy,
                                         y = y,
                                         muu = muu,
                                         mu = mu,
                                         transform = TRUE,
                                         obj_fun = elbo_fun,
                                         opt = opt,
                                         verbose = verbose)
        }

        if(xu_opt == "simultaneous")
        {
          results <- norm_grad_ascent_vi(cov_par_start = cov_par_start,
                                         cov_fun = "sqexp",
                                         dcov_fun_dtheta = dcov_fun_dtheta,
                                         dcov_fun_dknot = dcov_fun_dknot,
                                         knot_opt = 1:nrow(xu),
                                         xu = xu,
                                         xy = xy,
                                         y = y,
                                         muu = muu,
                                         mu = mu,
                                         transform = TRUE,
                                         obj_fun =  elbo_fun,
                                         opt = opt,
                                         verbose = verbose)
        }

        if(xu_opt == "oat")
        {
          results <- oat_knot_selection_norm_vi(cov_par_start = cov_par_start,
                                                cov_fun = cov_fun,
                                                dcov_fun_dtheta = dcov_fun_dtheta,
                                                dcov_fun_dknot = dcov_fun_dknot,
                                                xu_start = xu,
                                                xy = xy,
                                                y = y,
                                                mu = mu,
                                                muu_start = muu,
                                                obj_fun = elbo_fun,
                                                opt = opt,
                                                transform = TRUE,
                                                verbose = verbose,
                                                proposal_fun = knot_prop_ego_norm_vi,
                                                predict_laplace = predict_laplace
          )
        }

        if(xu_opt == "random")
        {
          results <- oat_knot_selection_norm_vi(cov_par_start = cov_par_start,
                                                cov_fun = cov_fun,
                                                dcov_fun_dtheta = dcov_fun_dtheta,
                                                dcov_fun_dknot = dcov_fun_dknot,
                                                xu_start = xu,
                                                xy = xy,
                                                y = y,
                                                mu = mu,
                                                muu_start = muu,
                                                obj_fun = elbo_fun,
                                                opt = opt,
                                                transform = TRUE,
                                                verbose = verbose,
                                                proposal_fun = knot_prop_random_norm_vi,
                                                predict_laplace = predict_laplace
          )
        }
      }

      if(vi == FALSE)
      {
        if(xu_opt == "fixed")
        {
          dcov_fun_dknot <- NA
          results <- norm_grad_ascent(cov_par_start = cov_par_start,
                                      cov_fun = cov_fun,
                                      dcov_fun_dtheta = dcov_fun_dtheta,
                                      dcov_fun_dknot = dcov_fun_dknot,
                                      knot_opt = NA, xu = xu,
                                      xy = xy,
                                      y = y,
                                      muu = muu,
                                      mu = mu,
                                      transform = TRUE,
                                      obj_fun = obj_fun_norm,
                                      opt = opt,
                                      verbose = verbose)
        }

        if(xu_opt == "simultaneous")
        {
          results <- norm_grad_ascent(cov_par_start = cov_par_start,
                                      cov_fun = "sqexp",
                                      dcov_fun_dtheta = dcov_fun_dtheta,
                                      dcov_fun_dknot = dcov_fun_dknot,
                                      knot_opt = 1:nrow(xu),
                                      xu = xu,
                                      xy = xy,
                                      y = y,
                                      muu = muu,
                                      mu = mu,
                                      transform = TRUE,
                                      obj_fun = obj_fun_norm,
                                      opt = opt,
                                      verbose = verbose)
        }

        if(xu_opt == "oat")
        {
          results <- oat_knot_selection_norm(cov_par_start = cov_par_start,
                                             cov_fun = cov_fun,
                                             dcov_fun_dtheta = dcov_fun_dtheta,
                                             dcov_fun_dknot = dcov_fun_dknot,
                                             xu_start = xu,
                                             xy = xy,
                                             y = y,
                                             mu = mu,
                                             muu_start = muu,
                                             obj_fun = obj_fun_norm,
                                             opt = opt,
                                             transform = TRUE,
                                             verbose = verbose,
                                             proposal_fun = knot_prop_ego_norm,
                                             predict_laplace = predict_laplace
          )
        }

        if(xu_opt == "random")
        {
          results <- oat_knot_selection_norm(cov_par_start = cov_par_start,
                                             cov_fun = cov_fun,
                                             dcov_fun_dtheta = dcov_fun_dtheta,
                                             dcov_fun_dknot = dcov_fun_dknot,
                                             xu_start = xu,
                                             xy = xy,
                                             y = y,
                                             mu = mu,
                                             muu_start = muu,
                                             obj_fun = obj_fun_norm,
                                             opt = opt,
                                             transform = TRUE,
                                             verbose = verbose,
                                             proposal_fun = knot_prop_random_norm,
                                             predict_laplace = predict_laplace
          )
        }
      }

    }
    ##########################
    ## Poisson case
    ##########################
    if(family == "poisson")
    {
      m <- args$m
      if(xu_opt == "fixed")
      {
        dcov_fun_dknot <- NA
        results <- laplace_grad_ascent(cov_par_start = cov_par_start,
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = NA,
                                       xu = xu,
                                       xy = xy,
                                       y = y,
                                       ff = rep(log(mean(y)), times = length(y)) - log(m),
                                       grad_loglik_fn = grad_loglik_fn_pois,
                                       dlog_py_dff = dlog_py_dff_pois,
                                       d2log_py_dff = d2log_py_dff_pois,
                                       d3log_py_dff = d3log_py_dff_pois,
                                       mu = mu,
                                       muu = muu,
                                       obj_fun = obj_fun_pois,
                                       opt = opt,
                                       verbose = verbose,
                                       predict_laplace = predict_laplace,
                                       m = m,
                                       ...)
      }
      if(xu_opt == "simultaneous")
      {
        results <- laplace_grad_ascent(cov_par_start = cov_par_start,
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = 1:nrow(xu),
                                       xu = xu,
                                       xy = xy,
                                       y = y,
                                       ff = rep(log(mean(y)), times = length(y)) - log(m),
                                       grad_loglik_fn = grad_loglik_fn_pois,
                                       dlog_py_dff = dlog_py_dff_pois,
                                       d2log_py_dff = d2log_py_dff_pois,
                                       d3log_py_dff = d3log_py_dff_pois,
                                       mu = mu,
                                       muu = muu,
                                       obj_fun = obj_fun_pois,
                                       opt = opt,
                                       verbose = verbose,
                                       predict_laplace = predict_laplace,
                                       m = m,
                                       ...)
      }
      if(xu_opt == "oat")
      {
        results <- oat_knot_selection(cov_par_start = cov_par_start,
                                      cov_fun = cov_fun,
                                      dcov_fun_dtheta = dcov_fun_dtheta,
                                      dcov_fun_dknot = dcov_fun_dknot,
                                      xu_start = xu,
                                      xy = xy,
                                      y = y,
                                      ff = rep(log(mean(y)), times = length(y)) - log(m),
                                      grad_loglik_fn = grad_loglik_fn_pois,
                                      dlog_py_dff = dlog_py_dff_pois,
                                      d2log_py_dff = d2log_py_dff_pois,
                                      d3log_py_dff = d3log_py_dff_pois,
                                      mu = mu,
                                      muu_start = muu,
                                      obj_fun = obj_fun_pois,
                                      opt = opt,
                                      verbose = verbose,
                                      proposal_fun = knot_prop_ego,
                                      predict_laplace = predict_laplace,
                                      m = m,
                                      ...)
      }
      if(xu_opt == "random")
      {
        results <- oat_knot_selection(cov_par_start = cov_par_start,
                                      cov_fun = cov_fun,
                                      dcov_fun_dtheta = dcov_fun_dtheta,
                                      dcov_fun_dknot = dcov_fun_dknot,
                                      xu_start = xu,
                                      xy = xy,
                                      y = y,
                                      ff = rep(log(mean(y)), times = length(y)) - log(m),
                                      grad_loglik_fn = grad_loglik_fn_pois,
                                      dlog_py_dff = dlog_py_dff_pois,
                                      d2log_py_dff = d2log_py_dff_pois,
                                      d3log_py_dff = d3log_py_dff_pois,
                                      mu = mu,
                                      muu_start = muu,
                                      obj_fun = obj_fun_pois,
                                      opt = opt,
                                      verbose = verbose,
                                      proposal_fun = knot_prop_random,
                                      predict_laplace = predict_laplace,
                                      m = m,
                                      ...)
      }
    }
    ##########################
    ## Bernoulli case
    ##########################
    if(family == "bernoulli")
    {
      if(xu_opt == "fixed")
      {
        dcov_fun_dknot <- NA
        results <- laplace_grad_ascent(cov_par_start = cov_par_start,
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = NA,
                                       xu = xu,
                                       xy = xy,
                                       y = y,
                                       ff = rep(0, times = length(y)),
                                       grad_loglik_fn = grad_loglik_fn_bern,
                                       dlog_py_dff = dlog_py_dff_bern,
                                       d2log_py_dff = d2log_py_dff_bern,
                                       d3log_py_dff = d3log_py_dff_bern,
                                       mu = mu,
                                       muu = muu,
                                       obj_fun = obj_fun_bern,
                                       opt = opt,
                                       verbose = verbose,
                                       predict_laplace = predict_laplace,
                                       ...)
      }
      if(xu_opt == "simultaneous")
      {
        results <- laplace_grad_ascent(cov_par_start = cov_par_start,
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = 1:nrow(xu),
                                       xu = xu,
                                       xy = xy,
                                       y = y,
                                       ff = rep(0, times = length(y)),
                                       grad_loglik_fn = grad_loglik_fn_bern,
                                       dlog_py_dff = dlog_py_dff_bern,
                                       d2log_py_dff = d2log_py_dff_bern,
                                       d3log_py_dff = d3log_py_dff_bern,
                                       mu = mu,
                                       muu = muu,
                                       obj_fun = obj_fun_bern,
                                       opt = opt,
                                       verbose = verbose,
                                       predict_laplace = predict_laplace,
                                       ...)
      }
      if(xu_opt == "oat")
      {
        results <- oat_knot_selection(cov_par_start = cov_par_start,
                                      cov_fun = cov_fun,
                                      dcov_fun_dtheta = dcov_fun_dtheta,
                                      dcov_fun_dknot = dcov_fun_dknot,
                                      xu_start = xu,
                                      xy = xy,
                                      y = y,
                                      ff = rep(0, times = length(y)),
                                      # ff = global_ff,
                                      grad_loglik_fn = grad_loglik_fn_bern,
                                      dlog_py_dff = dlog_py_dff_bern,
                                      d2log_py_dff = d2log_py_dff_bern,
                                      d3log_py_dff = d3log_py_dff_bern,
                                      mu = mu,
                                      muu_start = muu,
                                      obj_fun = obj_fun_bern,
                                      opt = opt,
                                      verbose = verbose,
                                      proposal_fun = knot_prop_ego,
                                      predict_laplace = predict_laplace,
                                      ...)
      }
      if(xu_opt == "random")
      {
        results <- oat_knot_selection(cov_par_start = cov_par_start,
                                      cov_fun = cov_fun,
                                      dcov_fun_dtheta = dcov_fun_dtheta,
                                      dcov_fun_dknot = dcov_fun_dknot,
                                      xu_start = xu,
                                      xy = xy,
                                      y = y,
                                      ff = rep(0, times = length(y)),
                                      # ff = global_ff,
                                      grad_loglik_fn = grad_loglik_fn_bern,
                                      dlog_py_dff = dlog_py_dff_bern,
                                      d2log_py_dff = d2log_py_dff_bern,
                                      d3log_py_dff = d3log_py_dff_bern,
                                      mu = mu,
                                      muu_start = muu,
                                      obj_fun = obj_fun_bern,
                                      opt = opt,
                                      verbose = verbose,
                                      proposal_fun = knot_prop_random,
                                      predict_laplace = predict_laplace,
                                      ...)
      }
    }
  }

  ## Full GP case
  if(sparse == FALSE)
  {
    if(family == "gaussian")
    {
      results <- norm_grad_ascent_full(cov_par_start = cov_par_start,
                                       cov_fun = cov_fun,
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       xy = xy,
                                       y = y,
                                       mu = mu,
                                       transform = TRUE,
                                       obj_fun = obj_fun_norm_full,
                                       opt = opt,
                                       verbose = verbose)
    }

    if(family == "poisson")
    {
      results <- laplace_grad_ascent_full(cov_par_start = cov_par_start,
                                          cov_fun = cov_fun,
                                          dcov_fun_dtheta = dcov_fun_dtheta,
                                          xy = xy,
                                          y = y,
                                          ff = rep(log(mean(y)), times = length(y)) - log(m),
                                          grad_loglik_fn = grad_loglik_fn_pois_full,
                                          dlog_py_dff = dlog_py_dff_pois,
                                          d2log_py_dff = d2log_py_dff_pois,
                                          d3log_py_dff = d3log_py_dff_pois,
                                          mu = mu,
                                          transform = TRUE,
                                          verbose = verbose,
                                          obj_fun = obj_fun_pois_full,
                                          opt = opt,
                                          m = m)
    }

    if(family == "bernoulli")
    {
      results <- laplace_grad_ascent_full(cov_par_start = cov_par_start,
                                          cov_fun = cov_fun,
                                          dcov_fun_dtheta = dcov_fun_dtheta,
                                          xy = xy,
                                          y = y,
                                          ff = rep(0, times = length(y)),
                                          grad_loglik_fn = grad_loglik_fn_bern_full,
                                          dlog_py_dff = dlog_py_dff_bern,
                                          d2log_py_dff = d2log_py_dff_bern,
                                          d3log_py_dff = d3log_py_dff_bern,
                                          mu = mu,
                                          transform = TRUE,
                                          verbose = verbose,
                                          obj_fun = obj_fun_bern_full,
                                          opt = opt)
    }
  }



  ## if the file path is a character string, save results there
  if(is.character(file_path))
  {
    saveRDS(object = list("results" = results,
                          "sparse" = sparse,
                          "family" = family,
                          "delta" = opt$delta,
                          "xu_init" = xu),
            file = file_path)
  }

  return(list("results" = results,
              "sparse" = sparse,
              "family" = family,
              "delta" = opt$delta,
              "xu_init" = xu))
}
