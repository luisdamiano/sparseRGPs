## GP VI functions required to implement GP VI approximation from
##  Titsias

## compute the trace term (as called that in Bauer et al. 2016)
#' Calculate the trace term in the ELBO.
#' @description Calculate the trace term in the ELBO for use in the variational approximation.
#'
#' @param cov_par A named list of covariance parameters.
#' @param Sigma12 Covariance matrix between observed data locations (rows) and knots (columns)
#' @param Sigma22 Variance covariance matrix at the knots.
#' @param delta A small diagonal perturbation to Sigma22 for numerical stability.
#' @return The value of the trace term in the ELBO.
#' @export
trace_term_fun <- function(cov_par, Sigma12, Sigma22, delta)
{
  tau <- cov_par$tau
  sigma <- cov_par$sigma
  Z2 <- solve(a = Sigma22, b = t(Sigma12))
  Z3 <- Sigma12 * t(Z2)
  Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)

  ## Note that Lambda here is the conditional variance matrix of the latent
  ##  function, so tau does not appear
  Lambda <- sigma^2 + delta - Z4

  return(-(1/(2 * tau^2)) * sum(Lambda))
}


## compute derivative of trace term
#' Derivative of the trace term.
#' @description Calculate the derivative of the trace term with respect to the nugget/noise.
#'
#' @param cov_par A named list of covariance parameters.
#' @param trace_term Scalar value of the current trace term in the ELBO.
#' @return The derivative of the trace term with respect to the nugget/noise.
#' @export
dtrace_term_dtau <- function(cov_par, trace_term)
{
  ## note that this is the derivative wrt log tau
  tau <- cov_par$tau
  sigma <- cov_par$sigma
  return(-2 * trace_term)
}

#' Derivative of the trace term.
#' @description Calculate the derivative of the trace term with respect to any covariance
#' parameters different from the nugget.
#'
#' @param cov_par A named list of covariance parameters.
#' @param A_trace A value calculated internally in the ELBO gradient function.
#' @return The derivative of the trace term with respect to the nugget/noise.
#' @export
dtrace_term_dcov_par <- function(cov_par, A_trace)
{
  tau <- cov_par$tau
  sigma <- cov_par$sigma

  return(-(1/(2 * tau^2)) * sum(A_trace))
}

## compute the elbo
#' @export
elbo_fun <- function(ff = NA, mu, Z, Sigma12,
                     Sigma22, y,
                     trace_term_fun = trace_term_fun,
                     cov_par, ...)
{
  ## ff are the (maximized) GP values at the observed data locations
  ## y is the vector of observed data
  ## m is a vector of the areas of each grid cell where counts are observed
  ## mu is a vector of GP means of the same length as ff
  ## Sigma12 is the matrix of covariances between the GP at the observed and knot locations
  ## Sigma22 is the variance covariance matrix at the knot locations
  ## Z is a vector of -----  sigma^2 + tau^2 - Sigma12 %*% Sigma22^(-1) %*% Sigma21
  ## trace_term_fun is a function to compute the trace term

  ## make sure ff is numeric
  # ff <- as.numeric(ff)
  y <- as.numeric(y)

  ## second derivative of p(y|xu, theta) wrt ff
  args <- list(...)
  delta <- args$delta

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

  ## compute the trace term
  trace_term <- trace_term_fun(cov_par = cov_par,
                               Sigma12 = Sigma12,
                               Sigma22 = Sigma22,
                               delta = delta)

  return(
    quad_form_part + det_part - (length(y) / 2) * log(2*pi) + trace_term
  )

}


## gradient of the ELBO wrt covariance parameters
#' @export
delbo_dcov_par <- function(cov_par,
                               cov_fun,
                               dcov_fun_dtheta,
                               dcov_fun_dknot = NA,
                               knot_opt,
                               xu, xy, y,
                               ff = NA,
                               mu, transform = TRUE,
                               delta = 1e-6, ...)
{
  ## cov_par is a list of the covariance function parameters
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is a string denoting either "sqexp" or "exp" for the squared exponential or exponential covariance functions
  ## xy are the observed data locations (matrix where rows are observations)
  ## xu are the unobserved knot locations (matrix where rows are knot locations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## ff is the current value for the ff vector that maximizes log(p(y|lambda)) + log(p(ff|theta, xu))
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
  ## dcov_fun_dknot is a single function that gives the derivative of the covariance function
  ##    wrt to the second input location x2, which I will always use for the knot location.
  ## knot_opt is a vector giving the knots over which to optimize, other knots are held constant
  ## transform is a logical argument that says whether to optimize on scales such that paramters are unconstrained
  ## ... are argument to be passed to the functions taking derivatives of log(p(y|lambda)) wrt ff

  ## store transformed parameter values if transform == TRUE
  if(transform == TRUE)
  {
    trans_par <- cov_par
  }

  ## initialize gradient
  if(is.list(dcov_fun_dtheta))
  {
    grad <- numeric(length = length(cov_par))
    names(grad) <- names(cov_par)
  }
  if(!is.list(dcov_fun_dtheta))
  {
    grad <- 0
  }

  if(is.function(dcov_fun_dknot))
  {
    grad_knot <- numeric(length = nrow(xu) * ncol(xu))

    ## get bounds for knots
    knot_bounds_l <- apply(X = xy, MARGIN = 2, FUN = min)
    knot_bounds_u <- apply(X = xy, MARGIN = 2, FUN = max)
    diffs <- knot_bounds_u - knot_bounds_l
    knot_bounds <- cbind(knot_bounds_l - diffs/10, knot_bounds_u + diffs/10)

  }
  trans_knot <- xu

  ## make sure y is a vector
  y <- as.numeric(y)

  ## create matrices that we will need for every derivative calculation
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")

    ## create Sigma12
    Sigma12 <- make_cov_mat_ardC(x = xy, x_pred = xu,
                                 cov_par = cov_par,
                                 cov_fun = cov_fun,
                                 delta = delta,
                                 lnames = lnames)

    ## create Sigma22
    Sigma22 <- make_cov_mat_ardC(x = xu,
                                 x_pred = matrix(),
                                 cov_fun = cov_fun,
                                 cov_par = cov_par,
                                 delta = delta,
                                 lnames = lnames) - as.list(cov_par)$tau^2 * diag(nrow(xu))
  }
  else{

    lnames <- character()

    ## create Sigma12
    Sigma12 <- make_cov_matC(x = xy, x_pred = xu,
                                 cov_par = cov_par,
                                 cov_fun = cov_fun,
                                 delta = delta)

    ## create Sigma22
    Sigma22 <- make_cov_matC(x = xu,
                                 x_pred = matrix(),
                                 cov_fun = cov_fun,
                                 cov_par = cov_par,
                                 delta = delta) - as.list(cov_par)$tau^2 * diag(nrow(xu))
  }

  ## create Z and W vectors and quantities
  ## create Z
  # Z2 <- solve(a = Sigma22, b = t(Sigma12))
  FF <- solve(a = Sigma22, b = t(Sigma12))
  Z <- numeric()
  Z <- rep(cov_par$tau^2 + delta, times = nrow(Sigma12))
  B <- 1/Z
  R <- chol(Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12))
  # FF <- solve(a = Sigma22, b = t(Sigma12))



  ## create the inverse of a commonly used m x m matrix (if m is the number of knot locations)
  ## maybe replace with cholesky decomposition eventually
  # RC <- chol(Sigma22 + t(Sigma12) %*% (B * Sigma12))
  C <- solve(a = Sigma22 + t(Sigma12) %*% (B * Sigma12))

  ## Perform other matrix computations shared by all derivatives

  ## compute component 2
  ## t(y - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (y - mu)
  comp2_1 <- (1/Z) * (y - mu) - t(solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*%
    (t((1/Z) * Sigma12) %*%  (y - mu))


  ## compute current trace term
  current_trace_term <- trace_term_fun(cov_par = cov_par,
                                   Sigma12 = Sigma12,
                                   Sigma22 = Sigma22,
                                   delta = delta)

  ## skip derivatives of the covariance function wrt theta if dcov_fun_dtheta is not a list
  if(is.list(dcov_fun_dtheta))
  {
    ## loop through each covariance parameter
    for(p in 1:length(cov_par))
    {
      par_name <- names(cov_par)[p]

      ## store transformed parameter values
      if(par_name %in% lnames)
      {
        trans_par[[p]] <- pos_to_real(x = cov_par[[which(names(cov_par) == par_name)]])
      }
      else{
        trans_par[[which(names(cov_par) == par_name)]] <-
          dcov_fun_dtheta[[which(names(dcov_fun_dtheta) == par_name)]](x1 = 0,
                                                                       x2 = 0,
                                                                       cov_par = cov_par,
                                                                       transform = transform)$trans_par

      }
      ## first compute dSigma12/dtheta_j
      if(cov_fun == "ard")
      {
        ## first compute dSigma12/dtheta_j
        dSigma12_dtheta <- dsig_dtheta_ardC(x = xy,
                                            x_pred = xu,
                                            cov_par = cov_par,
                                            cov_fun = cov_fun,
                                            par_name = par_name,
                                            lnames = lnames)

        ## next compute dSigma22/dtheta_j
        dSigma22_dtheta <- dsig_dtheta_ardC(x = xu,
                                            x_pred = matrix(),
                                            cov_par = cov_par,
                                            cov_fun = cov_fun,
                                            par_name = par_name,
                                            lnames = lnames)
      }

      else{
        ## first compute dSigma12/dtheta_j
        dSigma12_dtheta <- dsig_dthetaC(x = xy,
                                        x_pred = xu,
                                        cov_par = cov_par,
                                        cov_fun = cov_fun,
                                        par_name = par_name)

        ## next compute dSigma22/dtheta_j
        dSigma22_dtheta <- dsig_dthetaC(x = xu,
                                        x_pred = matrix(),
                                        cov_par = cov_par,
                                        cov_fun = cov_fun,
                                        par_name = par_name)
      }

      ## NOTE: for Gaussian data, tau is a special case where dSigma22_dtheta = 0
      if(par_name == "tau")
      {
        dSigma22_dtheta <- matrix(0, nrow = nrow(xu), ncol = nrow(xu))
      }

      ##########################################################################
      ## create components to compute the dervative of the trace term
      ## Create the vector A_trace = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
      ## Note that for the noise variance/tau, A is zero... At least in the Gaussian data case.
      if(par_name != "tau")
      {


        A1_trace <- numeric(length = nrow(xy))
        if(cov_fun != "ard" | !(par_name %in% lnames))
        {
          for(i in 1:length(A1_trace))
          {
            temp <- eval(
              substitute(expr = dcov_fun_dtheta$par(x1 = a, x2 = b, cov_par = c, transform = d),
                         env = list("a" = xy[i,], "b" = xy[i,], "c" = cov_par,
                                    "par" = par_name, "d" = transform))
            )
            A1_trace[i] <- temp$derivative
          }
        }


        #A1 <- rep(2 * cov_par$sigma * (ifelse(test = transform == TRUE, yes = cov_par$sigma, no = 1)), times = nrow(xy))
        A2_trace <- numeric(length = nrow(xy))
        # temp1 <- (2 * dSigma12_dtheta - Sigma12 %*% solve(a = Sigma22, b = dSigma22_dtheta))
        temp1 <- (2 * dSigma12_dtheta - t(FF) %*% dSigma22_dtheta )
        A2_trace <- apply(X = temp1 * t(FF), MARGIN = 1, FUN = sum)
        # temp2 <- FF
        # for(i in 1:length(A2))
        # {
        #   A2[i] <- temp1[i,] %*% temp2[,i]
        # }
        A_trace <- A1_trace - A2_trace
        ##########################################################################
      }


      ## Create the vector A = diag( dvar(y)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
      ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
      A1 <- numeric(length = nrow(xy))
      if(par_name == "tau")
      {
        for(i in 1:length(A1))
        {
          temp <- eval(
            substitute(expr = dcov_fun_dtheta$par(x1 = a, x2 = b, cov_par = c, transform = d),
                       env = list("a" = xy[i,], "b" = xy[i,], "c" = cov_par, "par" = par_name, "d" = transform))
          )
          A1[i] <- temp$derivative
        }
      }

      A <- A1

      ## we compute the gradient by summing some components together

      ## compute component 1
      ## trace( (Sigma)_inverse dSigma/dtheta_j)
      comp1_1 <- sum(A * B) - sum( diag( C %*% t(Sigma12) %*% (B * A * B * Sigma12)) )
      # comp1_1 <- sum(A * B) - sum( diag( solve( a = RC, b = solve( t(RC) , b = t(Sigma12) ) ) %*% (B * A * B * Sigma12)) )


      comp1_2_1 <- 2 * sum( diag(solve(a = Sigma22, b = t(Sigma12) %*% (B * dSigma12_dtheta))) )

      comp1_2_2 <- sum( diag( FF %*% t(solve(a = Sigma22, b = t(B * Sigma12))) %*% dSigma22_dtheta ))

      comp1_2_3 <- 2 * sum( diag( (FF %*% (B * Sigma12)) %*% (C %*% t(Sigma12) %*% (B * dSigma12_dtheta)) ))

      comp1_2_4 <- sum( diag( (FF %*% (B * Sigma12) ) %*%
                                ((C %*% t(Sigma12)) %*% (B * Sigma12) ) %*%
                                solve(a = Sigma22, b = dSigma22_dtheta)))

      comp1 <-  comp1_1 + comp1_2_1 - comp1_2_2 - comp1_2_3 + comp1_2_4

      ## compute component 2
      ## t(y - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (y - mu)
      # comp2_1 <- (1/Z) * (y - mu) - t(solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*%
      #   (t((1/Z) * Sigma12) %*%  (y - mu))
      comp2_2 <- A * comp2_1 + 2 * dSigma12_dtheta %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1) -
        t(FF) %*% dSigma22_dtheta %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1)

      comp2 <- t(comp2_1) %*% comp2_2

      ## compute derivative of the trace term
      if(par_name == "tau")
      {
        dtrace_term <- dtrace_term_dtau(cov_par = cov_par,
                                        trace_term = current_trace_term)
      }
      if(par_name != "tau")
      {
        dtrace_term <- dtrace_term_dcov_par(cov_par = cov_par,
                                            A_trace = A_trace)
      }

      ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
      ##    the p-th covariance parameter)
      grad[which(names(grad) == par_name)] <- (1/2) * comp2 -
        (1/2) * comp1 + dtrace_term

    }
  }


  ## optimize the knot locations if there is a derivative function for the knot locations
  ## loop through each covariance parameter
  if(is.function(dcov_fun_dknot))
  {
    ## make next two functions return sparse matrix
    ## function to create d(Sigma12)/dknot[k,d]
    dsig12_dknot <- function(k,d, cov_par,
                             dcov_fun_dknot,
                             xu, xy,
                             bounds = knot_bounds,
                             transform = TRUE)
    {

      ## m indexes the knot
      ## d indexes the dimension of the knot
      vec_m <- numeric(length = nrow(xy))
      for(i in 1:length(vec_m))
      {
        temp <- eval(
          substitute(expr = dcov_fun_dknot(x1 = a, x2 = b, cov_par = c, transform = e, bounds = h),
                     env = list("a" = xy[i,], "b" = xu[k,], "c" = cov_par, "e" = transform, "h" = bounds))
        )
        vec_m[i] <- temp$derivative[d]
      }

      mat <- Matrix::Matrix(data = 0, nrow = nrow(xy), ncol = nrow(xu))
      mat[,k] <- vec_m

      return(mat)
    }

    ## function to create d(Sigma22)/dknot[k,d]
    dsig22_dknot <- function(k,d,
                             cov_par,
                             dcov_fun_dknot,
                             xu, xy,
                             transform = TRUE,
                             bounds = knot_bounds)
    {
      ## note that because the derivative function is wrt x2, the fixed index m
      ##    must always be in the second argument
      mat <- Matrix::Matrix(data = 0, nrow = nrow(xu), ncol = nrow(xu))
      for(i in 1:nrow(mat))
      {
        temp <- eval(
          substitute(expr = dcov_fun_dknot(x1 = a, x2 = b, cov_par = c, transform = d, bounds = h),
                     env = list("a" = xu[i,], "b" = xu[k,], "c" = cov_par, "d" = transform, "h" = bounds))
        )
        mat[i,k] <- temp$derivative[d]
      }
      for(i in 1:nrow(mat))
      {
        temp <- eval(
          substitute(expr = dcov_fun_dknot(x1 = a, x2 = b, cov_par = c, transform = d, bounds = h),
                     env = list("a" = xu[i,], "b" = xu[k,], "c" = cov_par, "d" = transform, "h" = bounds))
        )
        mat[k,i] <- temp$derivative[d]
      }
      return(mat)
    }

    ## k indexes the knot, d indexes the dimension within the knot
    p <- 0
    # for(k in 1:length(knot_opt))
    for(k in 1:nrow(xu))
    {
      ## store transformed parameter values
      if(transform == TRUE)
      {
        trans_knot[k,] <- dcov_fun_dknot(x1 = 0, x2 = xu[k,], cov_par = cov_par,
                                         transform = transform, bounds = knot_bounds)$trans_par
      }
      for(d in 1:ncol(xu))
      {
        ## p is the index for the knot gradient vector of length m * d
        p <- p + 1

        ## make the gradient zero if current m is not one of the knots we
        ##  are optimizing over
        if(!k %in% knot_opt)
        {
          grad_knot[p] <- 0
        }
        if(k %in% knot_opt)
        {
          # par_name <- names(cov_par)[p]

          ## first compute dSigma12/dknot[k,d]
          dSigma12_dknot <- dsig12_dknot(k = k, d = d,
                                         cov_par = cov_par,
                                         dcov_fun_dknot = dcov_fun_dknot,
                                         xu = xu,
                                         xy = xy,
                                         transform = transform,
                                         bounds = knot_bounds)

          ## next compute dSigma22/dtheta_j
          dSigma22_dknot <- dsig22_dknot(k = k, d = d,
                                         cov_par = cov_par,
                                         dcov_fun_dknot = dcov_fun_dknot,
                                         xu = xu,
                                         transform = transform,
                                         bounds = knot_bounds)

          ##########################################################################
          ## create components to compute the dervative of the trace term
          ## Create the vector A_trace = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
          ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
          A1_trace <- 0

          #A1 <- rep(2 * cov_par$sigma * (ifelse(test = transform == TRUE, yes = cov_par$sigma, no = 1)), times = nrow(xy))
          A2_trace <- numeric(length = nrow(xy))
          temp1 <- (2 * dSigma12_dknot - t(FF) %*% dSigma22_dknot)
          A2_trace <- apply(X = temp1 * t(FF), MARGIN = 1, FUN = sum)
          # temp2 <- FF
          # for(i in 1:length(A2))
          # {
          #   A2[i] <- temp1[i,] %*% temp2[,i]
          # }
          A_trace <- A1_trace - A2_trace
          ##########################################################################


          ## Create the vector A = diag( dvar(y)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
          ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
          A1 <- 0
          A <- A1

          ## we compute the gradient by summing some components together

          ## compute component 1
          ## trace( (Sigma)_inverse dSigma/dtheta_j)
          comp1_1 <- sum(A * B) - sum( diag(C %*% t(Sigma12) %*% (B * A * B * Sigma12)) )

          comp1_2_1 <- 2 * sum( diag(solve(a = Sigma22, b = as.matrix(t(Sigma12) %*% (B * dSigma12_dknot))) ) )

          comp1_2_2 <- sum( diag( as.matrix(FF %*% t(solve(a = Sigma22, b = t(B * Sigma12))) %*% dSigma22_dknot) ))

          comp1_2_3 <- 2 * sum( diag( as.matrix( (FF %*% (B * Sigma12)) %*% (C %*% t(Sigma12) %*% (B * dSigma12_dknot))) ))

          comp1_2_4 <- sum( diag( (FF %*% (B * Sigma12) ) %*%
                                    ((C %*% t(Sigma12)) %*% (B * Sigma12) ) %*%
                                    solve(a = Sigma22, b = as.matrix(dSigma22_dknot) )))

          comp1 <-  comp1_1 + comp1_2_1 - comp1_2_2 - comp1_2_3 + comp1_2_4

          ## compute component 2
          ## t(y - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (y - mu)
          # comp2_1 <- (1/Z) * (y - mu) - t(solve(a = Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12), b = t( (1/Z) * Sigma12 ))) %*%
          #   (t((1/Z) * Sigma12) %*%  (y - mu))

          # comp2_1 <- (1/Z) * (y - mu) - t( solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*%
          #   (t((1/Z) * Sigma12) %*%  (y - mu))

          comp2_2 <- A * comp2_1 + 2 * dSigma12_dknot %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1) -
            t(FF) %*% dSigma22_dknot %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1)

          comp2 <- t(comp2_1) %*% comp2_2

          ## compute derivative of the trace term wrt the knot
          dtrace_term <- dtrace_term_dcov_par(cov_par = cov_par, A_trace = A_trace)

          ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
          ##    the p-th knot)
          grad_knot[p] <- (1/2) * as.numeric(comp2) -
            (1/2) * as.numeric(comp1) + dtrace_term
        }

      }
    }
  }
  if(is.function(dcov_fun_dknot))
  {
    return(list("gradient" = grad, "knot_gradient" = grad_knot, "trans_par" = trans_par, "trans_knot" = trans_knot))
  }
  if(is.na(dcov_fun_dknot))
  {
    return(list("gradient" = grad, "trans_par" = trans_par))
  }
}

## Gradient Ascent with the ELBO
#' @export
norm_grad_ascent_vi <- function(cov_par_start,
                                 cov_fun,
                                 dcov_fun_dtheta,
                                 dcov_fun_dknot,
                                 knot_opt,
                                 xu,
                                 xy,
                                 y,
                                 mu = NA,
                                 muu = NA,
                                 obj_fun = elbo_fun,
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

  ## This is an old argument that is necessary for gradient functions
  ##    indicating the need for the gradient of a transformed parameter
  transform <- TRUE

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
      xu[i,] <- xu[i,] + rnorm(n = ncol(xu),sd = 1e-5)
    }
  }

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
  optim_par = list("decay" = opt_master$decay, "epsilon" = opt_master$epsilon,
                   "eta" = opt_master$eta, "learn_rate" = opt_master$learn_rate)
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

  current_obj_fun <- obj_fun(mu = mu,
                             Z = Z,
                             Sigma12 = Sigma12,
                             Sigma22 = Sigma22,
                             y = y,
                             trace_term_fun = trace_term_fun,
                             cov_par = cov_par, delta = delta, ...) ## store the ELBO value

  obj_fun_vals[1] <- current_obj_fun
  cov_par_vals[1,] <- current_cov_par
  temp_grad_eval <- delbo_dcov_par(cov_par = cov_par, cov_fun = cov_fun, ## store the results of gradient function
                                       dcov_fun_dtheta = dcov_fun_dtheta,
                                       dcov_fun_dknot = dcov_fun_dknot,
                                       knot_opt = knot_opt,
                                       xu = xu, xy = xy,
                                       y = y, ff = NA, mu = mu, delta = delta,
                                       transform = transform)

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
          xu <- xu + optim_par$learn_rate * matrix(data = current_grad_val_knot,
                                                  nrow = nrow(xu),
                                                  ncol = ncol(xu),
                                                  byrow = TRUE)
          xu_vals <- abind::abind(xu_vals, xu, along = 3)
        }
      }

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
      Z <- numeric()
      Z <- rep(cov_par$tau^2 + delta, times = nrow(Sigma12))

      current_obj_fun <- obj_fun(mu = mu,
                                 Z = Z,
                                 Sigma12 = Sigma12,
                                 Sigma22 = Sigma22,
                                 y = y,
                                 trace_term_fun = trace_term_fun,
                                 cov_par = cov_par, delta = delta, ...) ## store the ELBO value

      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)

      temp_grad_eval <- delbo_dcov_par(cov_par = as.list(current_cov_par), cov_fun = cov_fun, ## store the results of gradient function
                                           dcov_fun_dtheta = dcov_fun_dtheta,
                                           dcov_fun_dknot = dcov_fun_dknot,
                                           knot_opt = knot_opt,
                                           xu = xu, xy = xy,
                                           y = y, ff = NA, mu = mu,
                                       delta = delta, transform = transform)
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
      # Z2 <- solve(a = Sigma22, b = t(Sigma12))
      Z <- numeric()
      Z <- rep(cov_par$tau^2 + delta, times = nrow(Sigma12))

      current_obj_fun <- obj_fun(mu = mu,
                                 Z = Z,
                                 Sigma12 = Sigma12,
                                 Sigma22 = Sigma22,
                                 y = y,
                                 trace_term_fun = trace_term_fun,
                                 cov_par = cov_par, delta = delta, ...) ## store the ELBO value
      obj_fun_vals[iter] <- current_obj_fun
      cov_par_vals <- rbind(cov_par_vals,current_cov_par)
      temp_grad_eval <- delbo_dcov_par(cov_par = as.list(current_cov_par), cov_fun = cov_fun, ## store the results of gradient function
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

## predict with the VI model
#' @export
predict_vi <- function(u_mean,
                        u_var,
                        xu,
                        x_pred,
                        cov_fun,
                        cov_par,
                        mu,
                        muu,
                        full_cov = FALSE,
                        family = "gaussian",
                        delta = 1e-6)
{
  ## u_mean is the posterior mean at the knot locations
  ## u_var is the posterior variance at the knot locations
  ## xu is a matrix where rows correspond to knots
  ## cov_fun is the covariance function
  ## cov_par are the covariance parameters
  ## mu is the marginal mean of the GP at locations x_pred
  ## full_cov is logical indicating whether full variance covariance matrix for predictions should be used

  ## create Sigma22
  ## create Sigma12 where the 1 corresponds to locations at which we'd like predictions
  if(family == "gaussian")
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(xu), sep = "")
      Sigma22 <- make_cov_mat_ardC(x = xu,
                               x_pred = matrix(),
                               cov_fun = cov_fun ,
                               cov_par = cov_par,
                               delta = delta,
                               lnames = lnames) -
        cov_par$tau^2 * diag(length(u_mean))
    }
    else{
      Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun , cov_par = cov_par, delta = delta) -
        cov_par$tau^2 * diag(length(u_mean))
    }

  }
  if(family != "gaussian")
  {
    return("Error: Only Gaussian data currently supported.")
    # Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun , cov_par = cov_par, delta = delta)
  }
  ## calculate Sigma22 inverse
  Sigma22_inv <- solve(a = Sigma22)


  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xu), sep = "")
    Sigma12 <- make_cov_mat_ardC(x = x_pred,
                                 x_pred = xu,
                                 cov_fun = cov_fun ,
                                 cov_par = cov_par,
                                 delta = delta,
                                 lnames = lnames)
  }
  else{
    Sigma12 <- make_cov_matC(x = x_pred, x_pred = xu, cov_fun = cov_fun , cov_par = cov_par, delta = delta)

  }

  ## calculate the predicted mean
  pred_mean <- mu + Sigma12 %*% solve(a = Sigma22, b = u_mean - muu)

  ## calculate the predictive variance
  ## create Sigma11 which is the prior variance covariance matrix at the
  ##    locations at which we desire predictions
  if(full_cov == TRUE)
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x_pred), sep = "")
      Sigma11 <- make_cov_mat_ardC(x = x_pred,
                               x_pred = matrix(),
                               cov_fun = cov_fun,
                               cov_par = cov_par,
                               delta = delta,
                               lnames = lnames)
    }
    else{
      Sigma11 <- make_cov_matC(x = x_pred,
                               x_pred = matrix(),
                               cov_fun = cov_fun,
                               cov_par = cov_par,
                               delta = delta)
    }

    pred_var <- Sigma11 + Sigma12 %*%
      (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv) %*%
      t(Sigma12)
  }
  if(full_cov == FALSE)
  {
    temp22 <- (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv)

    Sigma11 <- cov_par$tau^2 + cov_par$sigma^2 + delta
    pred_var <- Sigma11 + apply(X = Sigma12 *
      t((-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv) %*%
      t(Sigma12)), MARGIN = 1, FUN = sum )
    #
    # pred_var <- numeric()
    # for(i in 1:nrow(Sigma12))
    # {
    #   pred_var[i] <- cov_par$tau^2 + (t(Sigma12[i,]) %*% solve(a = Sigma22, b = t(t(Sigma12[i,]) ))) +
    #     t(Sigma12[i,]) %*% temp22 %*% t(t(Sigma12[i,]))
    # }
  }

  return(list("pred_mean" = pred_mean, "pred_var" = pred_var))

}

## OAT with the VI model
#' @export
oat_knot_selection_norm_vi <- function(cov_par_start,
                                    cov_fun,
                                    dcov_fun_dtheta,
                                    dcov_fun_dknot,
                                    obj_fun = elbo_fun,
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
                    "maxit" = 1000, "obj_tol" = 1e-3, "grad_tol" = 1,
                    "maxit_nr" = 1000, "delta" = 1e-6,
                    "tol_knot" = 1e-3, "maxknot" = 30, "TTmin" = 10, "TTmax" = 40,
                    "ego_cov_par" = list("sigma" = 3, "l" = 1, "tau" = 1e-3),
                    "ego_cov_fun" = "exp",
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
  opt <- norm_grad_ascent_vi(cov_par_start = cov_par_start,
                          cov_fun = cov_fun,
                          dcov_fun_dtheta = dcov_fun_dtheta,
                          dcov_fun_dknot = NA, xu = xu,
                          xy = xy, ff = NA, y = y, knot_opt = NA,
                          mu = mu, muu =  muu,
                          # transform = transform,
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
  while(numknots < maxknot && ifelse(test = length(obj_fun_vals) == 1, yes = TRUE,
                                     no = ifelse(test = current_obj_fun - obj_fun_vals[length(obj_fun_vals) - 1] > tol_knot,
                                                 yes = TRUE, no = FALSE)))
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
    opt <- norm_grad_ascent_vi(cov_par_start = cov_par,
                            cov_fun = cov_fun,
                            dcov_fun_dtheta = dcov_fun_dtheta,
                            dcov_fun_dknot = dcov_fun_dknot,
                            xu = rbind(xu, knot_prop),
                            knot_opt = (nrow(xu) + 1):(nrow(xu) + nrow(knot_prop)),
                            xy = xy, ff = NA, y = y,
                            mu = mu, muu = c(muu, muu[1]),
                            # transform = transform,
                            obj_fun = obj_fun,
                            opt = opt_master,
                            verbose = verbose, ...
    )
    ga_steps[iter] <- opt$iter
    current_obj_fun <- opt$obj_fun[length(opt$obj_fun)]
    obj_fun_vals[iter] <- opt$obj_fun[length(opt$obj_fun)]


    if(current_obj_fun - obj_fun_vals[length(obj_fun_vals) - 1] > tol_knot)
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

## variational inference knot proposal functions
## knot prop EGO
#' @export
knot_prop_ego_norm_vi <- function(norm_opt,
                               # obj_fun,
                               # grad_loglik_fn,
                               # d2log_py_dff,
                               # maxit_nr = 2000,
                               # tol_nr = 1e-5,
                               # y,
                               opt = list(), ...)
{
  ## norm_opt is a list of output returned from VI gradient ascent function
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
    # pseudo_xu <- rbind(xu, pseudo_prop[i,])
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
                               delta = delta, lnames = lnames) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }
    else{
      Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                               cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
      Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                               cov_par = norm_opt$cov_par,
                               cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

    }

    ## create Z
    # Z2 <- solve(a = Sigma22, b = t(Sigma12))
    # Z3 <- Sigma12 * t(Z2)
    # Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
    Z <- rep(norm_opt$cov_par$tau^2 + delta, times = nrow(norm_opt$xy))

    ## meta model y values
    obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z,
                                Sigma12 = Sigma12,
                                Sigma22 = Sigma22, y = y,
                                trace_term_fun = trace_term_fun,
                                cov_par = cov_par, delta = delta))
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
                                   cov_fun = norm_opt$cov_fun,
                                   delta = delta, lnames = lnames) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }
        else{
          Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
          Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }
        ## create Z
        # Z2 <- solve(a = Sigma22, b = t(Sigma12))
        # Z3 <- Sigma12 * t(Z2)
        # Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
        Z <-  rep(norm_opt$cov_par$tau^2 + delta, times = nrow(xy))

        ## meta model y values
        obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z,
                                    Sigma12 = Sigma12,
                                    Sigma22 = Sigma22,
                                    y = y,
                                    trace_term_fun = trace_term_fun,
                                    cov_par = cov_par,
                                    delta = delta))
      }
    }

    obj_fun_vals <- c(obj_fun_vals, obj_fun_eval)
    # obj_fun_vals <- c(obj_fun_eval)

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
                               cov_par = cov_par, delta = delta, lnames = lnames)

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
    # Z2 <- solve(a = Sigma22, b = t(Sigma12))
    # Z3 <- Sigma12 * t(Z2)
    # Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
    Z <- rep(norm_opt$cov_par$tau^2 + delta, times = nrow(norm_opt$xy))

    ## meta model y values
    obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z,
                                Sigma12 = Sigma12,
                                Sigma22 = Sigma22,
                                y = y,
                                trace_term_fun = trace_term_fun,
                                cov_par = cov_par,
                                delta = delta))

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
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun,
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
        # Z2 <- solve(a = Sigma22, b = t(Sigma12))
        # Z3 <- Sigma12 * t(Z2)
        # Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
        Z <- rep(norm_opt$cov_par$tau^2 + delta, times = nrow(norm_opt$xy))

        ## meta model y values
        obj_fun_eval <- try(obj_fun(mu = norm_opt$mu,
                                    Z = Z,
                                    Sigma12 = Sigma12,
                                    Sigma22 = Sigma22,
                                    y = y,
                                    trace_term_fun = trace_term_fun,
                                    cov_par = cov_par,
                                    delta = delta))
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

## Knot proposal based on a random subset of data locations
#' @export
knot_prop_random_norm_vi <- function(norm_opt,
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
        next
        # return()
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
                               cov_par = norm_opt$cov_par,
                               cov_fun = norm_opt$cov_fun,
                               delta = delta,
                               lnames = lnames)
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
    Z <- rep(norm_opt$cov_par$tau^2 + delta, times = nrow(norm_opt$xy))

    ## meta model y values
    obj_fun_eval <- try(obj_fun(mu = norm_opt$mu, Z = Z,
                                Sigma12 = Sigma12,
                                Sigma22 = Sigma22,
                                y = y, cov_par = cov_par,
                                trace_term_fun = trace_term_fun,
                                delta = delta))
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
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun,
                                   delta = delta, lnames = lnames)
          Sigma22 <- make_cov_mat_ardC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun,
                                   delta = delta, lnames = lnames) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))
        }
        else{
          Sigma12 <- make_cov_matC(x = norm_opt$xy, x_pred = pseudo_xu,
                                   cov_par = norm_opt$cov_par, cov_fun = norm_opt$cov_fun, delta = delta)
          Sigma22 <- make_cov_matC(x = pseudo_xu, x_pred = matrix(),
                                   cov_par = norm_opt$cov_par,
                                   cov_fun = norm_opt$cov_fun, delta = delta) - as.list(norm_opt$cov_par)$tau^2 * diag(nrow(pseudo_xu))

        }


        ## create Z
        # Z2 <- solve(a = Sigma22, b = t(Sigma12))
        # Z3 <- Sigma12 * t(Z2)
        # Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
        Z <- rep(norm_opt$cov_par$tau^2 + delta, times = nrow(norm_opt$xy))

        ## meta model y values
        obj_fun_eval <- try(obj_fun(mu = norm_opt$mu,
                                    Z = Z,
                                    Sigma12 = Sigma12,
                                    Sigma22 = Sigma22,
                                    y = y, cov_par = cov_par,
                                    trace_term_fun = trace_term_fun,
                                    delta = delta))
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

