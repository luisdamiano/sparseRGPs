## Functions to do prediction using the Laplace approximation 

predict_laplace <- function(u_mean, 
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
      Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(), 
                               cov_fun = cov_fun , 
                               cov_par = cov_par, 
                               delta = delta, lnames = lnames) - 
        cov_par$tau^2 * diag(length(u_mean))
    }
    else{
      Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), 
                               cov_fun = cov_fun , 
                               cov_par = cov_par, 
                               delta = delta) - 
        cov_par$tau^2 * diag(length(u_mean))
      
    }
    
  }
  if(family != "gaussian")
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(xu), sep = "")
      Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(), cov_fun = cov_fun ,
                               cov_par = cov_par, delta = delta, lnames = lnames)
      
    }
    else{
      Sigma22 <- make_cov_matC(x = xu,
                               x_pred = matrix(),
                               cov_fun = cov_fun ,
                               cov_par = cov_par, delta = delta)
      
    }
  }
  ## calculate Sigma22 inverse
  Sigma22_inv <- solve(a = Sigma22)
  
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xu), sep = "")
    Sigma12 <- make_cov_mat_ardC(x = x_pred, x_pred = xu, cov_fun = cov_fun , 
                                 cov_par = cov_par, 
                                 delta = delta, lnames = lnames)
    
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
                                   delta = delta, lnames = lnames)
      
    }
    else{
      Sigma11 <- make_cov_matC(x = x_pred, 
                               x_pred = matrix(), 
                               cov_fun = cov_fun, 
                               cov_par = cov_par, delta = delta)
      
    }
    Sigma11 <- diag(diag(Sigma11 - Sigma12 %*% solve(a = Sigma22, b = t(Sigma12)))) +
      (Sigma12 %*% solve(a = Sigma22, b = t(Sigma12)))
    pred_var <- Sigma11 + Sigma12 %*% 
                        (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv) %*% 
                        t(Sigma12)
  }
  if(full_cov == FALSE)
  {
    temp22 <- (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv)
    
    pred_var <- numeric()
    for(i in 1:nrow(Sigma12))
    {
      pred_var[i] <- cov_par$sigma^2 + cov_par$tau^2 + 
        t(Sigma12[i,]) %*% temp22 %*% t(t(Sigma12[i,]))
    }
  }

  return(list("pred_mean" = pred_mean, "pred_var" = pred_var))
  
}


## Full GP predictions from Laplace approximation
predict_laplace_full <- function(ff,
                                 xy,
                                 x_pred,
                                 L,
                                 W,
                                 cov_par,
                                 cov_fun,
                                 grad_log_py_ff,
                                 mu, 
                                 mu_pred, 
                                 full_cov = FALSE, 
                                 delta = 1e-6)
{
  ## ff is the approximate posterior mean from the Laplace approximation
  ## xy are the observed data locations
  ## x_pred are the locations at which we wish to predict
  ## W is the diagonal matrix of second derivatives of log(p(y|f)) wrt f
  ## L = chol(I + sqrt(-W) %*% Sigma %*% sqrt(-W))
  ## grad_log_py_ff is the gradient of log(p(y|f)) wrt f
  ## cov_fun is the covariance function
  ## cov_par are the covariance parameters 
  ## mu is the marginal mean of the GP at locations rbind(x_pred, xy)
  ## full_cov is logical indicating whether full variance covariance matrix for predictions should be used
  
  ## create Sigma12
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(x_pred), sep = "")
    Sigma12 <- make_cov_mat_ardC(x = x_pred, 
                                 x_pred = xy, 
                                 cov_fun = cov_fun , 
                                 cov_par = cov_par,
                                 delta = delta,
                                 lnames = lnames)
    
  }
  else{
    Sigma12 <- make_cov_matC(x = x_pred, x_pred = xy, cov_fun = cov_fun , cov_par = cov_par, delta = delta)
    
  }

  v <- solve(a = L, b = (sqrt(-diag(W)) %*% t(Sigma12)))
  
  ## calculate the predictive variance 
  ## create Sigma11 which is the prior variance covariance matrix at the locations at which we desire predictions
  if(full_cov == TRUE)
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x_pred), sep = "")
      Sigma11 <- make_cov_mat_ardC(x = x_pred, 
                                   x_pred = matrix(), 
                                   cov_fun = cov_fun, 
                                   cov_par = cov_par, 
                                   delta = delta, lnames = lnames)
    }
    else{
      Sigma11 <- make_cov_matC(x = x_pred, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)
      
    }
    pred_var <- Sigma11 - t(v) %*% v
  }
  if(full_cov == FALSE)
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x_pred), sep = "")
      Sigma11 <- make_cov_mat_ardC(x = x_pred, 
                                   x_pred = matrix(), 
                                   cov_fun = cov_fun, 
                                   cov_par = cov_par, 
                                   delta = delta, lnames = lnames)
      
    }
    else{
      Sigma11 <- make_cov_matC(x = x_pred, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)
      
    }
    pred_var <- diag(Sigma11 - t(v) %*% v)
    
    ## calculate the predicted mean
    pred_mean <- mu_pred + Sigma12 %*% grad_log_py_ff
  }
  
  return(list("pred_mean" = pred_mean, "pred_var" = pred_var))
  
}

# predict_laplace <- function(u_mean, 
#                             u_var, 
#                             xu, 
#                             x_pred, 
#                             cov_fun, 
#                             cov_par, 
#                             mu, 
#                             muu, 
#                             full_cov = FALSE, 
#                             family = "gaussian", 
#                             delta = 1e-6)
# {
#   ## u_mean is the posterior mean at the knot locations
#   ## u_var is the posterior variance at the knot locations
#   ## xu is a matrix where rows correspond to knots
#   ## cov_fun is the covariance function
#   ## cov_par are the covariance parameters 
#   ## mu is the marginal mean of the GP at locations rbind(x_pred, xu)
#   ## full_cov is logical indicating whether full variance covariance matrix for predictions should be used
#   
#   ## create Sigma22
#   ## create Sigma12 where the 1 corresponds to locations at which we'd like predictions
#   if(family == "gaussian")
#   {
#     Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun , cov_par = cov_par, delta = delta) - 
#       cov_par$tau^2 * diag(length(u_mean))
#     
#   }
#   if(family != "gaussian")
#   {
#     Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun , cov_par = cov_par, delta = delta)
#   }
#   ## calculate Sigma22 inverse
#   Sigma22_inv <- solve(a = Sigma22)
#   
#   Sigma12 <- make_cov_matC(x = x_pred, x_pred = xu, cov_fun = cov_fun , cov_par = cov_par, delta = delta)
#   
#   ## calculate the predicted mean
#   pred_mean <- mu + Sigma12 %*% solve(a = Sigma22, b = u_mean - muu)
#   
#   ## calculate the predictive variance 
#   ## create Sigma11 which is the prior variance covariance matrix at the locations at which we desire predictions
#   if(full_cov == TRUE)
#   {
#     Sigma11 <- make_cov_matC(x = x_pred, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)
#     pred_var <- Sigma11 + Sigma12 %*% (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv) %*% t(Sigma12)
#   }
#   if(full_cov == FALSE)
#   {
#     temp22 <- (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv)
#     
#     pred_var <- numeric()
#     for(i in 1:nrow(Sigma12))
#     {
#       pred_var[i] <- cov_par$sigma^2 + cov_par$tau^2 + 
#         t(Sigma12[i,]) %*% temp22 %*% t(t(Sigma12[i,]))
#     }
#   }
#   
#   return(list("pred_mean" = pred_mean, "pred_var" = pred_var))
#   
# }

## predictions for full GP with Gaussian data
predict_gp_full <- function(xy,
                            y,
                            x_pred, 
                            cov_fun, 
                            cov_par, 
                            mu, 
                            mu_pred,
                            full_cov = FALSE, 
                            delta = 1e-6)
{
  ## xy is a matrix where rows correspond to observed data input locations
  ## y is the observed data values
  ## cov_fun is the covariance function
  ## cov_par are the covariance parameters 
  ## mu is the marginal mean of the GP at data locations
  ## mu_pred is the marginal mean of the GP at x_pred
  ## full_cov is logical indicating whether full variance covariance matrix for predictions should be used
  
  ## check class of y
  y <- as.numeric(y)
  
  ## create covariance matrices
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(xy), sep = "")
    Sigma22 <- make_cov_mat_ardC(x = xy, x_pred = matrix(), 
                                 cov_par = cov_par, cov_fun = cov_fun, delta = delta, lnames = lnames)
    
    Sigma12 <- make_cov_mat_ardC(x = x_pred, x_pred = xy, cov_fun = cov_fun ,
                                 cov_par = cov_par, delta = delta, lnames = lnames)
  }
  else{
    Sigma22 <- make_cov_matC(x = xy, x_pred = matrix(), cov_par = cov_par, cov_fun = cov_fun, delta = delta)
    
    Sigma12 <- make_cov_matC(x = x_pred, x_pred = xy, cov_fun = cov_fun , cov_par = cov_par, delta = delta)
  }
  
  ## calculate the predicted mean
  pred_mean <- mu_pred + Sigma12 %*% solve(a = Sigma22, b = y - mu)
  
  ## calculate the predictive variance 
  ## create Sigma11 which is the prior variance covariance matrix at the locations at which we desire predictions
  if(full_cov == TRUE)
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x_pred), sep = "")
      Sigma11 <- make_cov_mat_ardC(x = x_pred, 
                                   x_pred = matrix(), 
                                   cov_fun = cov_fun,
                                   cov_par = cov_par,
                                   delta = delta, lnames = lnames)
      
    }
    else{
      Sigma11 <- make_cov_matC(x = x_pred, x_pred = matrix(), 
                               cov_fun = cov_fun, cov_par = cov_par, delta = delta)
      
    }
    pred_var <- Sigma11 - Sigma12 %*% solve(a = Sigma22, b = t(Sigma12))
  }
  if(full_cov == FALSE)
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(x_pred), sep = "")
      Sigma11 <- make_cov_mat_ardC(x = x_pred, 
                                   x_pred = matrix(), 
                                   cov_fun = cov_fun,
                                   cov_par = cov_par,
                                   delta = delta, lnames = lnames)
      
    }
    else{
      Sigma11 <- make_cov_matC(x = x_pred, x_pred = matrix(), 
                               cov_fun = cov_fun, cov_par = cov_par, delta = delta)
      
    }
    pred_var <- diag(Sigma11 - Sigma12 %*% solve(a = Sigma22, b = t(Sigma12)))
  }
  
  return(list("pred_mean" = pred_mean, "pred_var" = pred_var))
  
}

###################################################
## gp prediction wrapper function
###################################################
predict_gp <- function(mod, x_pred, mu_pred = NA, full_cov, vi = FALSE)
{
  ## mod is the object returned by a parameter estimation/knot selection function
  ## x_pred is the matrix of inputs at which we would like to predict
  ## full_cov is a logical argument dictating whether the full covariance matrix 
  ##    for predictions is to be returned
  
  
  
  family <- mod$family
  sparse <- mod$sparse
  m <- mod$results
  delta <- mod$delta
  
  ## print error if try to use VI with non gaussian data
  if(vi == TRUE && family != "gaussian")
  {
    return("Error: VI not supported for non-gaussian data.")
  }
  
  if(family == "poisson")
  {
    inv_link_fn <- function(x)
    {
      return(exp(x))
    }
  }
  if(family == "bernoulli")
  {
    inv_link_fn <- function(x)
    {
      return(1 / (1 + exp(-x)))
    }
  }
  if(family == "gaussian")
  {
    inv_link_fn <- NA
  }
  
  if(!is.matrix(x_pred))
  {
    print("Warning: x_pred must be a matrix. I'll try to make the conversion.")
    x_pred <- matrix(data = x_pred, ncol = 1)
  }
  
  if(any(is.na(mu_pred)))
  {
    print("Warnings: you did not define the mean of the GP at locations at which you wish to make predictions. Setting the mean to be zero.")
    mu_pred <- rep(0, times = nrow(x_pred))
  }
  
  ## use full gp predictions if results are from full GP
  if(sparse == FALSE)
  {
    ## gaussian family 
    if(family == "gaussian")
    {
      pred <- predict_gp_full(xy = m$xy, 
                              y = m$y, 
                              x_pred = x_pred, 
                              cov_fun = m$cov_fun, 
                              cov_par = m$cov_par, mu = m$mu, 
                              mu_pred = mu_pred, 
                              full_cov = full_cov, 
                              delta = delta)
    }
    
    ## poisson family
    if(family == "poisson")
    {
      pred <- predict_laplace_full(ff = m$fmax, 
                                   xy = m$xy, 
                                   x_pred = x_pred, 
                                   L = m$L, 
                                   W = m$W, 
                                   cov_par = m$cov_par, 
                                   cov_fun = m$cov_fun,
                                   grad_log_py_ff = m$grad_log_py_ff, 
                                   mu = m$mu, mu_pred = mu_pred, 
                                   full_cov = full_cov, 
                                   delta = delta)
    }
    
    ## bernoulli family
    if(family == "bernoulli")
    {
      pred <- predict_laplace_full(ff = m$fmax, 
                                   xy = m$xy, 
                                   x_pred = x_pred, 
                                   L = m$L, 
                                   W = m$W, 
                                   cov_par = m$cov_par, 
                                   cov_fun = m$cov_fun,
                                   grad_log_py_ff = m$grad_log_py_ff, 
                                   mu = m$mu, mu_pred = mu_pred, 
                                   full_cov = full_cov, 
                                   delta = delta)
    }
  }
  
  ## sparse predictions
  if(sparse == TRUE)
  {
    if(vi == TRUE)
    {
      pred <- predict_vi(u_mean = m$u_mean, 
                         u_var = m$u_var, 
                         xu = m$xu, 
                         x_pred = x_pred, 
                         cov_fun = m$cov_fun, 
                         cov_par = m$cov_par, 
                         mu = mu_pred, 
                         muu = m$muu, 
                         full_cov = full_cov, 
                         delta = delta)
    }
    if(vi == FALSE)
    {
      pred <- predict_laplace(u_mean = m$u_mean[1:nrow(m$xu)], 
                              u_var = m$u_var, 
                              xu = m$xu, 
                              x_pred = x_pred, 
                              cov_fun = m$cov_fun, 
                              cov_par = m$cov_par, 
                              mu = mu_pred, 
                              muu = m$muu[1:nrow(m$xu)], 
                              full_cov = full_cov, 
                              family = family, 
                              delta = delta)
    }

  }
  
  return(list("pred" = pred, "sparse" = sparse, "family" = family, "x_pred" = x_pred, "inverse_link" = inv_link_fn))
}

## subset of regressors / projected process prediction
predict_sor <- function(u_mean, 
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
      Sigma22 <- make_cov_mat_ardC(x = xu, x_pred = matrix(), 
                               cov_fun = cov_fun , 
                               cov_par = cov_par, 
                               delta = delta, lnames = lnames) - 
        cov_par$tau^2 * diag(length(u_mean))
    }
    else{
      Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun , cov_par = cov_par, delta = delta) - 
        cov_par$tau^2 * diag(length(u_mean))
    }
    
    
  }
  if(family != "gaussian")
  {
    if(cov_fun == "ard")
    {
      lnames <- paste("l", 1:ncol(xu), sep = "")
      Sigma22 <- make_cov_mat_ardC(x = xu, 
                                   x_pred = matrix(), 
                                   cov_fun = cov_fun , 
                                   cov_par = cov_par, 
                                   delta = delta, lnames = lnames)
      
    }
    else{
      Sigma22 <- make_cov_matC(x = xu, x_pred = matrix(), cov_fun = cov_fun , cov_par = cov_par, delta = delta)
      
    }
  }
  ## calculate Sigma22 inverse
  Sigma22_inv <- solve(a = Sigma22)
  
  if(cov_fun == "ard")
  {
    lnames <- paste("l", 1:ncol(x_pred), sep = "")
    Sigma12 <- make_cov_mat_ardC(x = x_pred, x_pred = xu, cov_fun = cov_fun , 
                                 cov_par = cov_par, delta = delta, lnames = lnames)
    
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
    # Sigma11 <- make_cov_matC(x = x_pred, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = delta)
    Sigma11 <- cov_par$tau^2 * diag(nrow(Sigma12)) + (Sigma12 %*% solve(a = Sigma22, b = t(Sigma12)))
    pred_var <- Sigma11 + Sigma12 %*% 
      (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv) %*% 
      t(Sigma12)
  }
  if(full_cov == FALSE)
  {
    temp22 <- (-Sigma22_inv + Sigma22_inv %*% u_var %*% Sigma22_inv)
    
    pred_var <- numeric()
    for(i in 1:nrow(Sigma12))
    {
      pred_var[i] <- cov_par$tau^2 + (t(Sigma12[i,]) %*% solve(a = Sigma22, b = t(t(Sigma12[i,]) ))) +
        t(Sigma12[i,]) %*% temp22 %*% t(t(Sigma12[i,]))
    }
  }
  
  return(list("pred_mean" = pred_mean, "pred_var" = pred_var))
  
}
