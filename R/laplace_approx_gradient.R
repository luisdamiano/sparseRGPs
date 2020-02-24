## functions to compute the gradient of log( q(y | theta, xu)) wrt theta and xu (the covariance parameters and inducing point locations)
## NOTE: package Matrix::Matrix is needed to encode sparse matrices and glmnet is needed to multiply them
# Rcpp::sourceCpp(file = "covariance_functions.cpp")

## NOTE: This gradient is calculated from 4 somewhat complicated components. The first two are identical to
##    the Gaussian data case (replacing f with y) and gradients check out splendidly in the Gaussian case
##    for all covariance parameters. Components 3 and 4 check out when comparing my computations which
##    avoid computing full n x n matrices, to those that do (meaning absolute differences in components are 
##    less than 10^-15). Yet, in the bernoulli case, covariance parameter gradients seem a little different from
##    those resulting from finite difference approximations. Note that I have checked the objective function
##    and newton raphson algorithms and they are fine. Also, the likelihood derivatives in the bernoulli case
##    are also correct. This means, either (a) one of the derivations is slightly incorrect for components 3 or 4
##    OR (b) there are numerical inaccuracies with all of the matrix multiplications which accrue. (a) might make sense 
##    if the gradient didn't consistently check out well in the Poisson case with sigma. The derivatives of covariance function wrt 
##    the length scale are also correct. Since the same exact code is used in both cases, I should be able to 
##    find a problem with the sigma gradients if there is actually a gradient problem, but I suppose that 
##    this should mean that numerical issues should be encountered in either case too (if that is indeed the issue).
##    It seems unlikely that the derivations are incorrect given how carefully I went over them
##    and that I carefully checked with others who did them as well. Also, the gradients wrt sigma in the Bernoulli data
##    case check out when the gradient is near zero, which is strange. All in all, everything is providing totally
##    reasonable and consistent results. It's hard to completely rule out the possibility of a small bug,
##    but I think that I am as close to that as humanly possible without carefully going over the derivations with someone
##    else willing.
dlogq_dcov_par <- function(cov_par, 
                           cov_fun, 
                           dcov_fun_dtheta,
                           dcov_fun_dknot = NA,
                           knot_opt,
                           xu, xy, y, 
                           ff, 
                           dlog_py_dff, 
                           d2log_py_dff, 
                           d3log_py_dff,
                           mu, transform = TRUE,
                           delta = 1e-6, ...)
{
  ## d2log_py_dff is a function that computes the vector of second derivatives 
  ##      of log(p(y|lambda)) wrt ff
  ## dlog_py_dff is a function that computes the vector of derivatives 
  ##      of log(p(y|lambda)) wrt ff
  ## cov_par is a list of the covariance function parameters 
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta unless ARD kernel is used
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

  
  ## create matrices that we will need for every derivative calculation
  ## change 1 for ard kernel
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
                                 lnames = lnames)
  }
  else{
    lnames <- character()
    
    
    ## create Sigma12
    Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
    
    ## create Sigma22
    Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)
  }
  
  
  ## create Z and W vectors and quantities
  # Z <- cov_par$sigma^2 + cov_par$tau^2 + delta - diag(Sigma12 %*% solve(a = Sigma22, b = t(Sigma12)))
  ## create Z
  # Z2 <- solve(a = Sigma22, b = t(Sigma12))
  FF <- solve(a = Sigma22, b = t(Sigma12))
  Z3 <- Sigma12 * t(FF)
  Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
  Z <- cov_par$sigma^2 + cov_par$tau^2 + delta - Z4
  
  W <- as.numeric(d2log_py_dff(ff = ff, y = y, ...))
  B <- 1/(Z - (1/W))
  R <- chol(Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12))
  # FF <- solve(a = Sigma22, b = t(Sigma12))
  W3 <- as.numeric(d3log_py_dff(ff = ff, y = y, ...))
  
  
  ## create the inverse of a commonly used m x m matrix (if m is the number of knot locations)
  ## maybe replace with cholesky decomposition eventually
  # RC <- chol(Sigma22 + t(Sigma12) %*% (B * Sigma12))
  C <- solve(a = Sigma22 + t(Sigma12) %*% (B * Sigma12))
  
  ## Perform other matrix computations shared by all derivatives 
  
  ## compute component 2
  ## t(ff - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (ff - mu)
  comp2_1 <- (1/Z) * (ff - mu) - t(solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*% 
    (t((1/Z) * Sigma12) %*%  (ff - mu))
  
  ## compute component 3 
  ## (I - Sigma %*% W)_inverse %*% dSigma/dtheta %*% dlogp(y|lambda)/dff
  grad_log_py_ff <- dlog_py_dff(y = y, ff = ff, ...)
  
  GG <- solve(a = Sigma22, b = t(Sigma12) %*% grad_log_py_ff)
  
  ## compute component 4 (see pages 18 and 28 of notes) this is actually a vector 
  ## THIS IS CORRECT
  ## diag( -W +Sigma_inverse)_inverse
  # E <- Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12) + t((1/Z) * Sigma12) %*% ((1/D) * (1/Z) * Sigma12)
  D <- W - 1/Z
  # E <- Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12) + t((1/Z) * Sigma12) %*% ((1/D) * (1/Z) * Sigma12)
  E2 <- diag(nrow(Sigma22)) + solve(a = t(R), b = t((1/Z) * Sigma12)) %*% t(solve(a = t(R), b = t(((1/D) * (1/Z) * Sigma12))))
  RE2 <- chol(E2)
  RE <- RE2 %*% R
  # RE <- chol(E)
  # Einv <- solve(a = E)
  
  comp4_1 <- numeric(length = nrow(xy))
  REinv <- solve(RE)
  for(i in 1:length(comp4_1))
  {
    temp_mat <- t(REinv) %*% (((1/Z) * (1/D))[i] * t(Sigma12)[,i]) 
    # temp_mat <- solve(a =  t(RE), b = (((1/Z) * (1/D))[i] * t(Sigma12)[,i]))
    # comp4_1[i] <- ((1/D)*(1/Z))[i] * Sigma12[i,] %*% Einv %*% (((1/Z) * (1/D))[i] * t(Sigma12)[,i])
    comp4_1[i] <- t(temp_mat) %*% temp_mat
  }
  comp4 <- -(1/D) + comp4_1
  
  
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
      
      
      
      ## Create the vector A = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
      ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
      A1 <- numeric(length = nrow(xy))
      
      if(cov_fun != "ard" | !(par_name %in% lnames))
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
      
      A2 <- numeric(length = nrow(xy))
      temp1 <- (2 * dSigma12_dtheta - t(FF) %*% dSigma22_dtheta )
      temp2 <- FF
      for(i in 1:length(A2))
      {
        A2[i] <- temp1[i,] %*% temp2[,i]
      }
      A <- A1 - A2
      
      ## we compute the gradient by summing some components together 
      
      ## compute component 1
      ## trace( (Sigma - W_inverse)_inverse dSigma/dtheta_j)
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
      ## t(ff - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (ff - mu)
      # comp2_1 <- (1/Z) * (ff - mu) - t(solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*% 
      #   (t((1/Z) * Sigma12) %*%  (ff - mu))
      comp2_2 <- A * comp2_1 + 2 * dSigma12_dtheta %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1) -
        t(FF) %*% dSigma22_dtheta %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1)
      
      comp2 <- t(comp2_1) %*% comp2_2
      
      ## compute component 3 
      ## (I - Sigma %*% W)_inverse %*% dSigma/dtheta %*% dlogp(y|lambda)/dff
      # grad_log_py_ff <- dlog_py_dff(y = y, ff, ...)
      
      ## comp3_1 looks right!!!
      ## This is also 4.1 in the notes
      comp3_1 <- A * grad_log_py_ff + 
        2 * dSigma12_dtheta %*% GG - 
        t(FF) %*% dSigma22_dtheta %*% GG
      
      ## This also checks out
      comp3 <- -(1/W) * (B * comp3_1) + 
        (1/W) * ((B * Sigma12) %*% (C %*% (t(Sigma12) %*% (B * comp3_1) )))
      
      
      # ## compute component 4 (see pages 18 and 28 of notes) this is actually a vector
      # ## diag( -W +Sigma_inverse)_inverse
      # # E <- Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12) + t((1/Z) * Sigma12) %*% ((1/D) * (1/Z) * Sigma12)
      # D <- W - 1/Z
      # # E <- Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12) + t((1/Z) * Sigma12) %*% ((1/D) * (1/Z) * Sigma12)
      # E2 <- diag(nrow(Sigma22)) + solve(a = t(R), b = t((1/Z) * Sigma12)) %*% t(solve(a = t(R), b = t(((1/D) * (1/Z) * Sigma12))))
      # RE2 <- chol(E2)
      # RE <- RE2 %*% R
      # # RE <- chol(E)
      # # Einv <- solve(a = E)
      # 
      # comp4_1 <- numeric(length = nrow(xy))
      # REinv <- solve(RE)
      # for(i in 1:length(comp4_1))
      # {
      #   temp_mat <- t(REinv) %*% (((1/Z) * (1/D))[i] * t(Sigma12)[,i]) 
      #   # temp_mat <- solve(a =  t(RE), b = (((1/Z) * (1/D))[i] * t(Sigma12)[,i]))
      #   # comp4_1[i] <- ((1/D)*(1/Z))[i] * Sigma12[i,] %*% Einv %*% (((1/Z) * (1/D))[i] * t(Sigma12)[,i])
      #   comp4_1[i] <- t(temp_mat) %*% temp_mat
      # }
      # comp4 <- -(1/D) + comp4_1
      
      ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
      ##    the p-th covariance parameter)
      grad[which(names(grad) == par_name)] <- (1/2) * comp2 -
        (1/2) * comp1 -
        (1/2) * (comp4 * (-W3)) %*% comp3
      
    }
  }
  
  
  ## optimize the knot locations if there is a derivative function for the knot locations
  ## loop through each covariance parameter
  if(is.function(dcov_fun_dknot))
  {
    ## make next two functions return sparse matrix
    ## function to create d(Sigma12)/dknot[k,d]
    dsig12_dknot <- function(k,d, 
                             cov_par, 
                             dcov_fun_dknot,
                             xu, 
                             xy, 
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
          
          ## Create the vector A = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
          ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
          A1 <- 0
          
          #A1 <- rep(2 * cov_par$sigma * (ifelse(test = transform == TRUE, yes = cov_par$sigma, no = 1)), times = nrow(xy))
          A2 <- numeric(length = nrow(xy))
          temp1 <- (2 * dSigma12_dknot - t(FF) %*% dSigma22_dknot)
          temp2 <- FF
          for(i in 1:length(A2))
          {
            A2[i] <- temp1[i,] %*% temp2[,i]
          }
          A <- A1 - A2
          
          ## we compute the gradient by summing some components together 
          
          ## compute component 1
          ## trace( (Sigma - W_inverse)_inverse dSigma/dtheta_j)
          comp1_1 <- sum(A * B) - sum( diag(C %*% t(Sigma12) %*% (B * A * B * Sigma12)) )
          
          ## The solve part of the following command weirdly throws an error 
          ## when b has class dgeMatrix
          comp1_2_1 <- 2 * sum( diag(solve(a = Sigma22, b = as.matrix(t(Sigma12) %*% (B * dSigma12_dknot))) ))
          
          ## This is now also behaving strangely
          comp1_2_2 <- sum( diag(as.matrix( FF %*% t(solve(a = Sigma22, b = t(B * Sigma12))) %*% dSigma22_dknot) ))
          
          comp1_2_3 <- 2 * sum( diag( as.matrix( (FF %*% (B * Sigma12)) %*% (C %*% t(Sigma12) %*% (B * dSigma12_dknot)) ) ))
          
          comp1_2_4 <- sum( diag( (FF %*% (B * Sigma12) ) %*% 
                                    ((C %*% t(Sigma12)) %*% (B * Sigma12) ) %*% 
                                    solve(a = Sigma22, b = as.matrix(dSigma22_dknot))))
          
          comp1 <-  comp1_1 + comp1_2_1 - comp1_2_2 - comp1_2_3 + comp1_2_4
          
          ## compute component 2
          ## t(ff - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (ff - mu)
          # comp2_1 <- (1/Z) * (ff - mu) - t(solve(a = Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12), b = t( (1/Z) * Sigma12 ))) %*% 
          #   (t((1/Z) * Sigma12) %*%  (ff - mu))
          
          # comp2_1 <- (1/Z) * (ff - mu) - t( solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*% 
          #   (t((1/Z) * Sigma12) %*%  (ff - mu))
          
          comp2_2 <- A * comp2_1 + 2 * dSigma12_dknot %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1) -
            t(FF) %*% dSigma22_dknot %*% solve(a = Sigma22, b = t(Sigma12) %*% comp2_1)
          
          comp2 <- t(comp2_1) %*% comp2_2
          
          ## compute component 3 
          ## (I - Sigma %*% W)_inverse %*% dSigma/dtheta %*% dlogp(y|lambda)/dff
          # grad_log_py_ff <- dlog_py_dff(y = y, ff, ...)
          
          comp3_1 <- A * grad_log_py_ff + 
            2 * dSigma12_dknot %*% GG - 
            t(FF) %*% dSigma22_dknot %*% GG
          
          comp3 <- -(1/W) * (B * comp3_1) + 
            (1/W) * ((B * Sigma12) %*% (C %*% (t(Sigma12) %*% (B * comp3_1) )))
          
          # ## compute component 4 (see pages 18 and 28 of notes) this is actually a vector
          # ## diag( -W +Sigma_inverse_inverse )_inverse
          # D <- W - 1/Z
          # # E <- Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12) + t((1/Z) * Sigma12) %*% ((1/D) * (1/Z) * Sigma12)
          # E2 <- diag(nrow(Sigma22)) + solve(a = t(R), b = t((1/Z) * Sigma12)) %*% t(solve(a = t(R), b = t(((1/D) * (1/Z) * Sigma12))))
          # RE2 <- chol(E2)
          # RE <- RE2 %*% R
          # # RE <- chol(E)
          # # Einv <- solve(a = E)
          # 
          # comp4_1 <- numeric(length = nrow(xy))
          # for(i in 1:length(comp4_1))
          # {
          #   # comp4_1[i] <- ((1/D)*(1/Z))[i] * Sigma12[i,] %*% Einv %*% (((1/Z) * (1/D))[i] * t(Sigma12)[,i])
          #   # comp4_1[i] <- t(solve(a =  t(RE), b = ((1/D)*(1/Z))[i] * Sigma12[i,])) %*% solve(a = t(RE), b = (((1/Z) * (1/D))[i] * t(Sigma12)[,i]))
          #   comp4_1[i] <- t(solve(a =  t(RE), b = (((1/Z) * (1/D))[i] * t(Sigma12)[,i]))) %*% solve(a = t(RE), b = (((1/Z) * (1/D))[i] * t(Sigma12)[,i]))
          #   
          # }
          # comp4 <- -(1/D) + comp4_1
          
          ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
          ##    the p-th knot)
          # print(class(comp1))
          # print(class(comp2))
          # print(class(comp3))
          # print(class(comp4))
          # print(class(W3))
          
          
          grad_knot[p] <- (1/2) * as.numeric(comp2) -
            (1/2) * as.numeric(comp1) -
            (1/2) * (as.numeric(comp4) * (-W3)) %*% as.numeric(comp3)
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


## FULL GP WITH LAPLACE APPROXIMATION
dlogq_dcov_par_full <- function(cov_par, 
                           cov_fun, 
                           dcov_fun_dtheta,
                           xy, y, 
                           ff, 
                           dlog_py_dff, 
                           d2log_py_dff, 
                           d3log_py_dff,
                           mu, transform = TRUE,
                           delta = 1e-6, ...)
{
  ## d3log_py_dff is a function that computes the vector of third derivatives 
  ##      of log(p(y|lambda)) wrt ff
  ## d2log_py_dff is a function that computes the vector of second derivatives 
  ##      of log(p(y|lambda)) wrt ff
  ## dlog_py_dff is a function that computes the vector of derivatives 
  ##      of log(p(y|lambda)) wrt ff
  ## cov_par is a list of the covariance function parameters 
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is a string denoting either "sqexp" or "exp" for the squared exponential or exponential covariance functions
  ## xy are the observed data locations (matrix where rows are observations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## ff is the current value for the ff vector that maximizes log(p(y|lambda)) + log(p(ff|theta, xu))
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in 
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
  ## transform is a logical argument that says whether to optimize on scales such that paramters are unconstrained
  ## ... are argument to be passed to the functions taking derivatives of log(p(y|lambda)) wrt ff 
  
  ## store transformed parameter values if transform == TRUE
  if(transform == TRUE)
  {
    trans_par <- cov_par
  }
  
  ## initialize gradient 
  grad <- numeric(length = length(cov_par))
  names(grad) <- names(cov_par)
  
  
  ## create matrices that we will need for every derivative calculation
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
    Sigma11 <- make_cov_mat_ardC(x = xy, x_pred = matrix(), 
                             cov_fun = cov_fun, 
                             cov_par = cov_par, 
                             delta = delta,
                             lnames = lnames)
  }
  else{
    lnames <- character()
    Sigma11 <- make_cov_matC(x = xy, x_pred = matrix(), 
                             cov_fun = cov_fun, 
                             cov_par = cov_par, 
                             delta = delta)
  }

  W <- as.numeric(d2log_py_dff(ff = ff, y = y, ...))
  W3 <- as.numeric(d3log_py_dff(ff = ff, y = y, ...))
  
  ## compute cholesky decomposition of I + sqrt(-W) %*% Sigma11 %*% sqrt(-W)
  L <- t(chol(diag(length(y)) + sqrt(-diag(W)) %*% Sigma11 %*% sqrt(-diag(W))))
  
  ## compute sqrt(-W) %*% L^T_inverse (L_inverse sqrt(-W))
  R <- sqrt(-(W)) * solve(a = t(L),
                                b = solve(a = L, b = sqrt(-diag(W))))
  
  ## compute L_inverse %*% (sqrt(-W) %*% Sigma11)
  C <- solve(a = L, b = (sqrt(-diag(W)) %*% Sigma11))
  
  ## compute -1/2 diag( diag(Sigma11) - diag(C^T C)) %*% d3log_py_dff
  s2 <- -1/2 * diag(
    diag(Sigma11)  - diag(t(C) %*% C)
  ) %*% (-W3)
  
  grad_log_py_ff <- dlog_py_dff(y = y, ff = ff, ...)
  
  ## calculate (-W) %*% f + (d/df) log(p(y|f))
  b <- (-diag(W)) %*% (ff - mu) + grad_log_py_ff
  
  ## calculate b - ( sqrt(-W) %*% t(L) )_inverse %*% (L_inverse sqrt(-W) %*% Sigma b)
  a <- b - sqrt(-diag(W)) %*% solve(a = t(L) , 
                                    b = solve(a = L, b = sqrt(-diag(W)) %*% Sigma11 %*% b))
  
  ## skip derivatives of the covariance function wrt theta if dcov_fun_dtheta is not a list
  ## loop through each covariance parameter
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
  }
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
    
    ## next compute dSigma11/dtheta_j
    # dSigma11_dtheta <- dsig22_dtheta(par_name = par_name, 
    #                                  cov_par = cov_par, 
    #                                  dcov_fun_dtheta = dcov_fun_dtheta, 
    #                                  xu = xu, 
    #                                  transform = transform)
    if(cov_fun == "ard")
    {
      dSigma11_dtheta <- dsig_dtheta_ardC(x = xy, 
                                      x_pred = matrix(), 
                                      cov_par = cov_par, 
                                      cov_fun = cov_fun, 
                                      par_name = par_name,
                                      lnames = lnames)
    }
    else{
      dSigma11_dtheta <- dsig_dthetaC(x = xy, 
                                      x_pred = matrix(), 
                                      cov_par = cov_par, 
                                      cov_fun = cov_fun, 
                                      par_name = par_name)
    }
    
    
    
    ## compute 1/2 a^T %*% dSigma11_dtheta %*% a - 1/2 tr(R dSigma11_dtheta)
    s1 <- 1/2 * (t(a) %*% dSigma11_dtheta %*% a) - 1/2 * sum(diag(R %*% dSigma11_dtheta))
    # s1 <- 1/2 * (t(solve(a = Sigma11, b = ff)) %*% dSigma11_dtheta %*% solve(a = Sigma11, b = ff)) - 1/2 * sum(diag(R %*% dSigma11_dtheta))
    
    ## compute dSigma11_dtheta %*% dlog_py_dff
    b <- dSigma11_dtheta %*% grad_log_py_ff
    
    ## compute b - Sigma11 %*% R %*% b
    s3 <- b - Sigma11 %*% R %*% b
    
    
    ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
    ##    the p-th covariance parameter)
    grad[which(names(grad) == par_name)] <- s1 + t(s2) %*% s3
    
  }
  
  return(list("gradient" = grad, "trans_par" = trans_par))
}

## TRY USING THE PRODUCT DENSITY FROM DETERIMANENTAL POINT PROCESS AS KNOT PRIOR

## functions to compute the gradient of log( p(y | theta, xu)) wrt theta and xu (the covariance parameters and inducing point locations)
## NOTE: package Matrix::Matrix is needed to encode sparse matrices and glmnet is needed to multiply them
## NOTE: This function is for Gaussian data
# Rcpp::sourceCpp(file = "covariance_functions.cpp")
dlogp_dcov_par <- function(cov_par, 
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
    Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
    
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
  Z3 <- Sigma12 * t(FF)
  Z4 <- apply(X = Z3, MARGIN = 1, FUN = sum)
  Z <- cov_par$sigma^2 + cov_par$tau^2 + delta - Z4
  
  B <- 1/Z
  R <- chol(Sigma22 + t(Sigma12) %*% ((1/Z) * Sigma12))
  
  
  
  ## create the inverse of a commonly used m x m matrix (if m is the number of knot locations)
  ## maybe replace with cholesky decomposition eventually
  # RC <- chol(Sigma22 + t(Sigma12) %*% (B * Sigma12))
  C <- solve(a = Sigma22 + t(Sigma12) %*% (B * Sigma12))
  
  ## Perform other matrix computations shared by all derivatives 
  
  ## compute component 2
  ## t(y - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (y - mu)
  comp2_1 <- (1/Z) * (y - mu) - t(solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*% 
    (t((1/Z) * Sigma12) %*%  (y - mu))
  
  
  ## skip derivatives of the covariance function wrt theta if dcov_fun_dtheta is not a list
  if(is.list(dcov_fun_dtheta))
  {
    ## loop through each covariance parameter
    for(p in 1:length(cov_par))
    {
      par_name <- names(cov_par)[p]
      
      # ## store transformed parameter values
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
      
      ## Create the vector A = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
      ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
      A1 <- numeric(length = nrow(xy))
      if(cov_fun != "ard" | !(par_name %in% lnames))
      {
        for(i in 1:length(A1))
        {
          temp <- eval(
            substitute(expr = dcov_fun_dtheta$par(x1 = a, x2 = b, cov_par = c, transform = d),
                       env = list("a" = xy[i,], "b" = xy[i,], "c" = cov_par, 
                                  "par" = par_name, "d" = transform))
          )
          A1[i] <- temp$derivative
        }
      }
      
      #A1 <- rep(2 * cov_par$sigma * (ifelse(test = transform == TRUE, yes = cov_par$sigma, no = 1)), times = nrow(xy))
      A2 <- numeric(length = nrow(xy))
      # temp1 <- (2 * dSigma12_dtheta - Sigma12 %*% solve(a = Sigma22, b = dSigma22_dtheta))
      temp1 <- (2 * dSigma12_dtheta - t(FF) %*% dSigma22_dtheta )
      temp2 <- FF
      for(i in 1:length(A2))
      {
        A2[i] <- temp1[i,] %*% temp2[,i]
      }
      A <- A1 - A2
      
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
      
      ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
      ##    the p-th covariance parameter)
      grad[which(names(grad) == par_name)] <- (1/2) * comp2 -
        (1/2) * comp1
      
    }
  }
  
  
  ## optimize the knot locations if there is a derivative function for the knot locations
  ## loop through each covariance parameter
  if(is.function(dcov_fun_dknot))
  {
    ## make next two functions return sparse matrix
    ## function to create d(Sigma12)/dknot[k,d]
    dsig12_dknot <- function(k,d, cov_par, dcov_fun_dknot, 
                             xu, xy, bounds = knot_bounds, 
                             transform = TRUE)
    {
      
      ## k indexes the knot
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
    dsig22_dknot <- function(k,d, cov_par, 
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
          
          ## Create the vector A = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
          ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
          A1 <- 0
          
          #A1 <- rep(2 * cov_par$sigma * (ifelse(test = transform == TRUE, yes = cov_par$sigma, no = 1)), times = nrow(xy))
          A2 <- numeric(length = nrow(xy))
          temp1 <- (2 * dSigma12_dknot - t(FF) %*% dSigma22_dknot)
          temp2 <- FF
          for(i in 1:length(A2))
          {
            A2[i] <- temp1[i,] %*% temp2[,i]
          }
          A <- A1 - A2
          
          ## we compute the gradient by summing some components together 
          
          ## compute component 1
          ## trace( (Sigma)_inverse dSigma/dtheta_j)
          comp1_1 <- sum(A * B) - sum( diag(C %*% t(Sigma12) %*% (B * A * B * Sigma12)) )
          
          comp1_2_1 <- 2 * sum( diag(solve(a = Sigma22, b = as.matrix(t(Sigma12) %*% (B * dSigma12_dknot)) )) )
          
          comp1_2_2 <- sum( diag( as.matrix(FF %*% t(solve(a = Sigma22, b = t(B * Sigma12))) %*% dSigma22_dknot) ))
          
          comp1_2_3 <- 2 * sum( diag( as.matrix((FF %*% (B * Sigma12)) %*% (C %*% t(Sigma12) %*% (B * dSigma12_dknot))) ))
          
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
          
          
          ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
          ##    the p-th knot)
          grad_knot[p] <- (1/2) * as.numeric(comp2) -
            (1/2) * as.numeric(comp1)
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

## NOTE: This function is fora full GP with Gaussian data
# Rcpp::sourceCpp(file = "covariance_functions.cpp")
dlogp_dcov_par_full <- function(cov_par, 
                           cov_fun, 
                           dcov_fun_dtheta,
                           xy, y, 
                           mu, transform = TRUE,
                           delta = 1e-6, ...)
{
  ## cov_par is a list of the covariance function parameters 
  ##    NAMES AND ORDERING MUST MATCH dcov_fun_dtheta
  ## cov_fun is a string denoting either "sqexp" or "exp" for the squared exponential or exponential covariance functions
  ## xy are the observed data locations (matrix where rows are observations)
  ## y is the vector of the observed data values
  ## mu is the mean of the GP at each of the observed data locations
  ## dcov_fun_dtheta is a list of functions which corresponds (in order) to the covariance parameters in 
  ##    cov_par. These functions compute the derivatives of the covariance function wrt the particular covariance parameter.
  ##    NAMES AND ORDERING MUST MATCH cov_par
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
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
  }
  
  ## make sure y is a vector
  y <- as.numeric(y)
  
  ## create matrices that we will need for every derivative calculation
  ## create Sigma11
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
    Sigma11 <- make_cov_mat_ardC(x = xy, 
                                        x_pred = matrix(), 
                                        cov_par = cov_par, 
                                        cov_fun = cov_fun, 
                                        lnames = lnames,
                                 delta = delta)
  }
  else{
    lnames <- character()
    Sigma11 <- make_cov_matC(x = xy, 
                                 x_pred = matrix(), 
                                 cov_par = cov_par, 
                                 cov_fun = cov_fun, 
                                 delta = delta)
  }

  ## Perform other matrix computations shared by all derivatives 
  
  ## alpha = (Sigma_11)_inverse %*% y 
  alpha <- solve(a = Sigma11, b = y)
  
  ## skip derivatives of the covariance function wrt theta if dcov_fun_dtheta is not a list
  ## loop through each covariance parameter
  for(p in 1:length(cov_par))
  {
    par_name <- names(cov_par)[p]
    
    # ## store transformed parameter values
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
    
    

    
    ## compute dSigma11/dtheta_j
    ## NOTE: for Gaussian data, tau is a special case where dSigma11_dtheta = 0
    if(cov_fun == "ard")
    {
      dSigma11_dtheta <- dsig_dtheta_ardC(x = xy, 
                                          x_pred = matrix(), 
                                          cov_par = cov_par, 
                                          cov_fun = cov_fun, 
                                          par_name = par_name,
                                          lnames = lnames)
    }
    else{
      dSigma11_dtheta <- dsig_dthetaC(x = xy, 
                                      x_pred = matrix(), 
                                      cov_par = cov_par, 
                                      cov_fun = cov_fun, 
                                      par_name = par_name)
    }
    
    
    
    
    
    ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
    ##    the p-th covariance parameter)
    grad[which(names(grad) == par_name)] <- (1/2) * sum(
      diag(
        (alpha %*% t(alpha)) %*% dSigma11_dtheta - solve(a = Sigma11, b = dSigma11_dtheta) 
      )
    )
    
  }
  
  return(list("gradient" = grad, "trans_par" = trans_par))
  
}



# ## test the gradient function
# y <- c(rpois(n = 2, lambda = c(2,2)))
# cov_par <- list("sigma" = 1.25, "l" = 2.27, "tau" = 1e-5)
# cov_fun <- mv_cov_fun_sqrd_exp
# dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl)
# dcov_fun_dknot <- dsqexp_dx2
# xy <- matrix(rep(c(1,2), times = 1), nrow = 2)
# xu <- matrix(rep(c(1), times = 1), nrow = 1)
# ff <- log(y)
# mu <- c(0,0)
# m <- c(1,1)
# 
## looks like this is probably working... currently takes 5.26 seconds for n = 1000, k = 121
# system.time(test_grad <- dlogq_dcov_par(cov_par = cov_par,
#                             cov_fun = "sqexp",
#                             dcov_fun_dtheta = dcov_fun_dtheta,
#                             dcov_fun_dknot = NA,
#                             xu = xu,
#                             xy = xy,
#                             y = y,
#                             ff = fstart,
#                             dlog_py_dff = dlog_py_dff_pois,
#                             d2log_py_dff = d2log_py_dff_pois,
#                             mu = mu,
#                             muu = muu,
#                             transform = TRUE,
#                             m = m))
# test_grad

## test making the derivative matrix
# dsig12_dtheta <- function(par_name, cov_par, dcov_fun_dtheta, xu, xy, transform = TRUE)
# {
#   mat <- matrix(nrow = nrow(xy), ncol = nrow(xu))
#   for(i in 1:nrow(mat))
#   {
#     for(j in 1:ncol(mat))
#     {
#       mat[i,j] <- eval(
#         substitute(expr = dcov_fun_dtheta$par(x1 = a, x2 = b, cov_par = c, transform = d),
#                    env = list("a" = xy[i,], "b" = xu[j,], "c" = cov_par, "par" = par_name, "d" = transform))
#       )
#     }
#   }
#   return(mat)
# }
# xy <- matrix(rep(c(1,2), times = 1), nrow = 2)
# xu <- matrix(rep(c(1.1,1.6), times = 1), nrow = 2)
# dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl)
# 
# dsig12_dtheta(par_name = "sigma", cov_par = list("sigma" = 1, "l" = 1), dcov_fun_dtheta = dcov_fun_dtheta, xu = xu, xy = xy,transform = TRUE)
# dsqexp_dsigma(x1 = xy[1,], x2 = xu[2,], cov_par = list("sigma" = 1, "l" = 1), transform = TRUE)
# 
# dsig12_dtheta(par_name = "l", cov_par = list("sigma" = 1, "l" = 1), dcov_fun_dtheta = dcov_fun_dtheta, xu = xu, xy = xy,transform = TRUE)
# dsqexp_dl(x1 = xy[1,], x2 = xu[2,], cov_par = list("sigma" = 1, "l" = 1), transform = TRUE)


## Function to take in covariance parameter values and knot locations and output the derivative of the laplace approximation...
  ## needs to interface with optim() so that BFGS optimization can be used. 
  ## This means that this function must do the Newton-Raphson optimization of the latent function
gradient_logq <- function(cov_par, ...)
{
  ## cov_par here must be a vector, not a list! But the elements of the vector should be obtained by 
    ## c(unlist(cov_par), as.numeric(t(xu))) for the cov_par list as passed to the other optimization functions
  ## ... optional arguments should include 
  ## cov_fun, 
  ## dcov_fun_dtheta,
  ## dcov_fun_dknot = NA,
  ## xu 
  ## xy 
  ## y 
  ## ff, 
  ## dlog_py_dff, 
  ## d2log_py_dff, 
  ## mu, 
  ## muu
  ## tol_nr
  ## maxit
  ## transform
  ## grad_loglik_fn
  ## obj_fun
  ## dlogq_dcov_par the function to compute gradients of the laplace approximation as used in the gradient ascent function
  ## num_cov_par is the number of elements in cov_par that correspond to covariance parameters that are not knots
  ##  these parameters must always come first in cov_par
  ## as well as other optional arguments that will be passed to the derivative functions of the data density
  
  args <- list(...)
  # cov_fun <- args$cov_fun
  # dcov_fun_dtheta <- args$dcov_fun_dtheta
  # xu <- args$xu
  # xy <- args$xy
  # y <- args$y
  # ff <- args$ff
  # dlog_py_dff <- args$dlog_py_dff
  # d2log_py_dff <- args$d2log_py_dff
  # mu <- args$mu
  # muu <- args$muu
  # transform = args$transform
  # grad_loglik_fn <- args$grad_loglik_fn
  # tol_nr <- args$tol_nr
  # maxit <- args$maxit
  # num_cov_par <- args$num_cov_par
  # dlogq_dcov_par <- args$dlogq_dcov_par
  
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
                              # obj_fun = args$obj_fun, 
                              # grad_loglik_fn = args$grad_loglik_fn, 
                              # d2log_py_dff = args$d2log_py_dff,
                              # tol = args$tol_nr,
                              # maxit = args$maxit,
                              cov_par = as.list(c(trans_cov_par[1:args$num_cov_par], args$other_cov_par)),
                              # cov_fun = args$cov_fun,
                              # xy = args$xy, 
                              xu = xu, 
                              # y = args$y, mu = args$mu, muu = args$muu, 
                              ...)
  
  ## compute the gradient 
  temp_grad_eval <- dlogq_dcov_par(cov_par = as.list(c(trans_cov_par[1:args$num_cov_par], args$other_cov_par)), ## store the results of the gradient function
                                    # cov_fun = args$cov_fun,
                                    # dcov_fun_dtheta = args$dcov_fun_dtheta,
                                    # dcov_fun_dknot = args$dcov_fun_dknot,
                                    xu = xu,
                                    # xy = args$xy,
                                    # y = args$y,
                                    ff = nr_step$gp,
                                    # dlog_py_dff = args$dlog_py_dff,
                                    # d2log_py_dff = args$d2log_py_dff,
                                    # mu = args$mu,
                                    # transform = args$transform,
                                   ...)
  
  # print(c(trans_cov_par[1:args$num_cov_par], other_cov_par))
  # print(c(temp_grad_eval$grad, temp_grad_eval$knot_gradient))
  
  if(is.function(dcov_fun_dknot))
  {
    grad <- temp_grad_eval$grad
    grad_knot <- temp_grad_eval$knot_gradient
    return(c(grad, grad_knot))
  }
  if(is.na(dcov_fun_dknot))
  {
    grad <- temp_grad_eval$grad
    return(grad)
  }
  
}

## SOR likelihood gradient
dlogp_dcov_par_sor <- function(cov_par, 
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
  ## create Sigma12
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
    Sigma12 <- make_cov_matC(x = xy, x_pred = xu, cov_par = cov_par, cov_fun = cov_fun, delta = delta)
    
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
  FF <- solve(a = Sigma22, b = t(Sigma12))
  
  
  
  ## create the inverse of a commonly used m x m matrix (if m is the number of knot locations)
  ## maybe replace with cholesky decomposition eventually
  # RC <- chol(Sigma22 + t(Sigma12) %*% (B * Sigma12))
  C <- solve(a = Sigma22 + t(Sigma12) %*% (B * Sigma12))
  
  ## Perform other matrix computations shared by all derivatives 
  
  ## compute component 2
  ## t(y - mu) %*% Sigma_inverse %*% dSigma/dtheta %*% Sigma_inverse %*% (y - mu)
  comp2_1 <- (1/Z) * (y - mu) - t(solve(a = R, b = solve(a = t(R), b = t( (1/Z) * Sigma12 )))) %*% 
    (t((1/Z) * Sigma12) %*%  (y - mu))
  
  
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
        lnames <- character()
        
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
      
      ## Create the vector A = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
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
      
      ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
      ##    the p-th covariance parameter)
      grad[which(names(grad) == par_name)] <- (1/2) * comp2 -
        (1/2) * comp1
      
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
          
          ## Create the vector A = diag( dvar(f)/dtheta - (d/dtheta) Sigma12 %*% Sigma22_inverse %*% Sigma21)
          ## if we are dealing with transformed parameters make sure you have the derivative dsigma/dsigma_transformed!!
          A1 <- 0
          A <- A1
          
          ## we compute the gradient by summing some components together 
          
          ## compute component 1
          ## trace( (Sigma)_inverse dSigma/dtheta_j)
          comp1_1 <- sum(A * B) - sum( diag(C %*% t(Sigma12) %*% (B * A * B * Sigma12)) )
          
          comp1_2_1 <- 2 * sum( diag(solve(a = Sigma22, b = as.matrix(t(Sigma12) %*% (B * dSigma12_dknot)) )) )
          
          comp1_2_2 <- sum( diag( as.matrix(FF %*% t(solve(a = Sigma22, b = t(B * Sigma12))) %*% dSigma22_dknot) ))
          
          comp1_2_3 <- 2 * sum( diag( as.matrix((FF %*% (B * Sigma12)) %*% (C %*% t(Sigma12) %*% (B * dSigma12_dknot))) ))
          
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
          
          
          ## compute dlog(q(y|theta))/dtheta (the derivative of the laplace approximation wrt
          ##    the p-th knot)
          grad_knot[p] <- (1/2) * as.numeric(comp2) -
            (1/2) * as.numeric(comp1)
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