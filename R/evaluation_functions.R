# Evaluation and assessment functions

## rmse function to true GP
my_rmse <- function(pred, actual)
{
  return(sqrt(mean((pred - actual)^2)))
}

## negative log probability function
my_nlp <- function(y, pred_mean, pred_var, family = "gaussian", par, mc = FALSE, mv = FALSE)
{
  if(family == "gaussian")
  {
    tau <- par$tau
    if(mv == FALSE)
    {
      nlp <- -dnorm(x = y, mean = pred_mean, sd = sqrt(pred_var + tau^2), log = TRUE)
    }
    if(mv == TRUE)
    {
      nlp <- -mvtnorm::dmvnorm(x = y, mean = pred_mean, 
                               sigma = pred_var + tau^2 * diag(nrow(pred_mean)),
                               log = TRUE)
    }
  }
  if(family == "bernoulli")
  {
    g <- function(ff, ...)
    {
      args <- list(...)
      y <- args$y
      pred_mean <- args$pred_mean
      pred_var <- args$pred_var
      p <- my_logistic(ff)
      ll <- dbinom(x = y, size = 1, prob = p, log = FALSE)
      return(
        ll * dnorm(x = ff, mean = pred_mean, sd = sqrt(pred_var), log = FALSE)
      )
    }
    if(mv == FALSE & mc == FALSE)
    {
      nlp <- numeric()
      
      for(i in 1:length(y))
      {
          nlp[i] <- -log(integrate(f = g, lower = -Inf, upper = Inf, 
                            y = y[i], pred_mean = as.numeric(pred_mean)[i],
                            pred_var = pred_var[i], par = par)$value)
      } 
    }
    if(mv == FALSE & mc == TRUE)
    {
      
      nlp <- numeric()
      mcsd <- numeric()
      for(i in 1:length(y))
      {
        ffsim <- rnorm(n = par$n, mean = pred_mean[i], sd = sqrt(pred_var[i]))
        ll <- dbinom(x = y[i], size = 1, prob = my_logistic(x = ffsim), log = FALSE)
        nlp[i] <- -log(mean(ll))
        mcsd[i] <- sd(ll)
      } 
    }
    if(mv == TRUE)
    {
      ## this does a monte carlo estimate
      ffsim <- mvtnorm::rmvnorm(n = par$n,
                                mean = pred_mean, sigma = pred_var)
      ymat <- matrix(rep(x = y, times = par$n), ncol = length(y), byrow = TRUE)
      mat <- abind::abind(ffsim, ymat, along = 3)
      g <- function(par_vec)
      {
        exp(sum(dbinom(x = par_vec[,2], size = 1, prob = my_logistic(par_vec[,1]), log = TRUE)))
      }
      ll <- apply(X = mat, MARGIN = 1, FUN = g)
      nlp <- -log(mean(ll))
      mcsd <- sd(ll)
    }
  }
  if(family == "poisson")
  {
    g <- function(ff, ...)
    {
      args <- list(...)
      y <- args$y
      pred_mean <- args$pred_mean
      pred_var <- args$pred_var
      par <- args$par
      m <- par$m
      ll <- dpois(x = y, lambda = exp(ff) * m, log = FALSE)
      return(
        ll * dnorm(x = ff, mean = pred_mean, sd = sqrt(pred_var), log = FALSE)
      )
    }
    
    if(mc == FALSE & mv == FALSE)
    {
      # nlp <- -log(integrate(f = g, lower = -Inf, upper = Inf,
      #                       y = y, pred_mean = pred_mean, pred_var = pred_var, par = par)$value)
      nlp <- numeric()
      
      for(i in 1:length(y))
      {
        nlp[i] <- -log(integrate(f = g, lower = -Inf, upper = Inf, 
                                 y = y[i], pred_mean = as.numeric(pred_mean)[i],
                                 pred_var = pred_var[i], par = par)$value)
      } 
    }
    
    if(mc == TRUE & mv == FALSE)
    {
      return("Error: Monte Carlo estimate of marginal NLP for Poisson data must be fixed.")
      # ffsim <- rnorm(n = par$n, mean = pred_mean, sd = sqrt(pred_var))
      # ll <- dpois(x = y, lambda = exp(ffsim) * par$m, log = FALSE)
      # nlp <- -log(mean(ll))
      # mcsd <- sd(ll)
    }
    if(mv == TRUE)
    {
      ffsim <- mvtnorm::rmvnorm(n = par$n,
                                mean = pred_mean, sigma = pred_var)
      ymat <- matrix(rep(x = y, times = par$n), ncol = length(y), byrow = TRUE)
      mat <- abind::abind(ffsim, ymat, along = 3)
      g <- function(par_vec, m)
      {
        exp(sum(dpois(x = par_vec[,2], size = 1, prob = exp(par_vec[,1]) * m, log = TRUE)))
      }
      ll <- apply(X = mat, MARGIN = 1, FUN = g, m = par$m)
      nlp <- -log(mean(ll))
      mcsd <- sd(ll)
    }
  }
  
  if(mc == FALSE)
  {
    return(list("nlp" = nlp))
  }
  if(mc == TRUE)
  {
    return(list("nlp" = nlp, "mcsd" = mcsd, "mcn" = par$n))
  }

}


## function calculate the KL divergence between two Multivariate Gaussians wrt 
##    the first density
my_kl <- function(mean1, mean2, sigma1, sigma2, mv = FALSE)
{
  if(mv == TRUE)
  {
    chol1 <- t(chol(x = sigma1))
    chol2 <- t(chol(x = sigma2))
    return(
      (1/2) * (
        sum(diag(solve(a = sigma2, b = sigma1))) + t(mean2 - mean1) %*% 
          solve(a = sigma2, b = mean2 - mean1) - length(mean1) + 
          sum(log(diag(chol2))) - sum(log(diag(chol1)))
      )
    )
  }
  if(mv == FALSE)
  {
    return(
      (1/2) * (
        sum(sigma1 / sigma2) + t(mean2 - mean1) %*% 
          ( diag((1/sigma2 )) %*% (mean2 - mean1)) - length(mean1) + 
          sum(log(sigma2)) - sum(log(sigma1))
      )
    )
  }
  
}

## function to take data, and results from predict_laplace, and produce a graph
pred_plot_fun <- function(y, xy, xu, xu_init, pred_results, ci_level, title = "Predictions", 
                          family = "gaussian", alpha = 0.2, sparse = FALSE, size = 0.2, 
                          xlabel = "x", ylabel = "y", samples = FALSE, ...)
{
  if(samples == FALSE)
  {
    quant <- qnorm(p = 1 - (1 - ci_level)/2, mean = 0, sd = 1)
    lb <- pred_results$pred_mean - quant * sqrt(pred_results$pred_var)
    ub <- pred_results$pred_mean + quant * sqrt(pred_results$pred_var)
    
    if(family == "gaussian")
    {
      pred_mean <- pred_results$pred_mean
    }
    if(family == "bernoulli")
    {
      pred_mean <- my_logistic(pred_results$pred_mean)
    }
    if(family == "poisson")
    {
      pred_mean <- exp(pred_results$pred_mean)
    }
  }
  if(samples == TRUE)
  {
    lb <- apply(X = pred_results$pred_samples, MARGIN = 2, 
                FUN = quantile, probs = (1 - ci_level)/2)
    ub <- apply(X = pred_results$pred_samples, MARGIN = 2, 
                FUN = quantile, probs = 1 - (1 - ci_level)/2)
    
    if(family == "gaussian")
    {
      pred_mean <- apply(X = pred_results$pred_samples, MARGIN = 2, 
                         FUN = mean)
    }
    if(family == "bernoulli")
    {
      pred_mean <- apply(X = my_logistic(pred_results$pred_samples), MARGIN = 2, 
                         FUN = mean)
    }
    if(family == "poisson")
    {
      pred_mean <- apply(X = exp(pred_results$pred_samples), MARGIN = 2, 
                         FUN = mean)
    }
    
  }
  
  if(sparse == TRUE)
  {
    if(family == "gaussian")
    {
      quant <- qnorm(p = 1 - (1 - ci_level)/2, mean = 0, sd = 1)
      mymin <- min(c(y, lb), na.rm = TRUE)
      mymax <- max(c(y, ub), na.rm = TRUE)
      small_int_y <- 2 * (mymax - mymin) / 100
      small_int_x <- (max(xy) - min(xy)) / 100
      
      myplot <- ggplot() + 
        geom_point(mapping = aes(x = xy, y = y), size = size) +
        geom_line(mapping = aes(x = xy, y = pred_mean), colour = "blue", linetype = "dashed") +
        geom_ribbon(mapping = aes(x = xy, 
                                  ymin = lb,
                                  ymax = ub),
                    fill = "blue", alpha = alpha) +
        # geom_line(mapping = aes(x = xy, y = lb), colour = "blue", linetype = 2) +
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        geom_segment(mapping = aes(x = xu, xend = xu,
                                   y = mymin - small_int_y, yend = mymin + small_int_y), col = "blue") +
        geom_segment(mapping = aes(x = xu - small_int_x, xend = xu + small_int_x,
                                   y = mymin, yend = mymin) , col = "blue") +
        geom_segment(mapping = aes(x = xu_init, xend = xu_init,
                                   y = mymax - small_int_y, yend = mymax + small_int_y), col = "red") +
        geom_segment(mapping = aes(x = xu_init - small_int_x, xend = xu_init + small_int_x,
                                   y = mymax, yend = mymax) , col = "red") +
        ggtitle(title)
    }
    if(family == "poisson")
    {
      args <- list(...)
      m <- args$m
      
      quant <- qnorm(p = 1 - (1 - ci_level)/2, mean = 0, sd = 1)
      mymin <- min(c(y, exp(lb) * m ), na.rm = TRUE)
      mymax <- max(c(y, exp(ub) * m ), na.rm = TRUE)
      small_int_y <- 2 * (mymax - mymin) / 100
      small_int_x <- (max(xy) - min(xy)) / 100
      
      myplot <- ggplot() + 
        geom_point(mapping = aes(x = xy, y = y), size = size) +
        geom_line(mapping = aes(x = xy, y = pred_mean * m), colour = "blue", linetype = "dashed") +
        geom_ribbon(mapping = aes(x = xy, 
                                  ymin = exp(lb) * m,
                                  ymax = exp(ub) * m), fill = "blue", alpha = alpha) +
        # geom_line(mapping = aes(x = xy, y = exp(lb) * m ), colour = "blue", linetype = 2) +
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        geom_segment(mapping = aes(x = xu, xend = xu,
                                   y = mymin - 3 * small_int_y, yend = mymin - small_int_y), col = "blue") +
        geom_segment(mapping = aes(x = xu - small_int_x, xend = xu + small_int_x,
                                   y = mymin - 2 * small_int_y, yend = mymin - 2 * small_int_y) , col = "blue") +
        geom_segment(mapping = aes(x = xu_init, xend = xu_init,
                                   y = mymax - small_int_y, yend = mymax + small_int_y), col = "red") +
        geom_segment(mapping = aes(x = xu_init - small_int_x, xend = xu_init + small_int_x,
                                   y = mymax, yend = mymax) , col = "red") +
        ggtitle(title)
    }
    
    if(family == "bernoulli")
    {
      my_logistic <- function(x)
      {
        return(1 / (1 + exp(-x)))
      }
      args <- list(...)
      
      quant <- qnorm(p = 1 - (1 - ci_level)/2, mean = 0, sd = 1)
      mymin <- min(c(y, my_logistic(lb)), na.rm = TRUE)
      mymax <- max(c(y, my_logistic(ub)), na.rm = TRUE) + 0.025
      small_int_y <- 2 * (mymax - mymin) / 100
      small_int_x <- (max(xy) - min(xy)) / 100
      
      myplot <- ggplot() + 
        geom_point(mapping = aes(x = xy, y = y), size = size) +
        geom_line(mapping = aes(x = xy, y = pred_mean), colour = "blue", linetype = "dashed") +
        geom_ribbon(mapping = aes(x = xy, 
                                  ymin = my_logistic(lb),
                                  ymax = my_logistic(ub)), fill = "blue", alpha = alpha) +
        # geom_line(mapping = aes(x = xy, y = my_logistic(lb) * m ), colour = "blue", linetype = 2) +
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        geom_segment(mapping = aes(x = xu, xend = xu,
                                   y = mymin - 3 * small_int_y, yend = mymin - small_int_y), col = "blue") +
        geom_segment(mapping = aes(x = xu - small_int_x, xend = xu + small_int_x,
                                   y = mymin - 2 * small_int_y, yend = mymin - 2 * small_int_y) , col = "blue") +
        geom_segment(mapping = aes(x = xu_init, xend = xu_init,
                                   y = mymax - small_int_y, yend = mymax + small_int_y), col = "red") +
        geom_segment(mapping = aes(x = xu_init - small_int_x, xend = xu_init + small_int_x,
                                   y = mymax, yend = mymax) , col = "red") +
        ggtitle(title)
    }
  }
  
  if(sparse == FALSE)
  {
    if(family == "gaussian")
    {
      quant <- qnorm(p = 1 - (1 - ci_level)/2, mean = 0, sd = 1)
      mymin <- min(c(y, lb), na.rm = TRUE)
      mymax <- max(c(y, ub), na.rm = TRUE)
      small_int_y <- 2 * ((mymax - mymin) / 100)
      small_int_x <- (max(xy) - min(xy)) / 100
      
      myplot <- ggplot() + 
        geom_point(mapping = aes(x = xy, y = y), size = size) +
        geom_line(mapping = aes(x = xy, y = pred_mean), colour = "blue", linetype = "dashed") +
        geom_ribbon(mapping = aes(x = xy, 
                                  ymin = lb,
                                  ymax = ub),
                    fill = "blue", alpha = alpha) +
        # geom_line(mapping = aes(x = xy, y = lb), colour = "blue", linetype = 2) +
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        ggtitle(title)
    }
    if(family == "poisson")
    {
      args <- list(...)
      m <- args$m
      
      quant <- qnorm(p = 1 - (1 - ci_level)/2, mean = 0, sd = 1)
      mymin <- min(c(y, exp(lb) * m ), na.rm = TRUE)
      mymax <- max(c(y, exp(ub) * m ), na.rm = TRUE)
      small_int_y <- 2 * (mymax - mymin) / 100
      small_int_x <- (max(xy) - min(xy)) / 100
      
      myplot <- ggplot() + 
        geom_point(mapping = aes(x = xy, y = y), size = size) +
        geom_line(mapping = aes(x = xy, y = pred_mean * m), colour = "blue", linetype = "dashed") +
        geom_ribbon(mapping = aes(x = xy, 
                                  ymin = exp(lb) * m,
                                  ymax = exp(ub) * m), fill = "blue", alpha = alpha) +
        # geom_line(mapping = aes(x = xy, y = exp(lb) * m ), colour = "blue", linetype = 2) +
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        ggtitle(title)
    }
    
    if(family == "bernoulli")
    {
      my_logistic <- function(x)
      {
        return(1 / (1 + exp(-x)))
      }
      args <- list(...)
      
      quant <- qnorm(p = 1 - (1 - ci_level)/2, mean = 0, sd = 1)
      mymin <- min(c(y, my_logistic(lb)), na.rm = TRUE)
      mymax <- max(c(y, my_logistic(ub)), na.rm = TRUE)
      small_int_y <- 2 * (mymax - mymin) / 100
      small_int_x <- (max(xy) - min(xy)) / 100
      
      myplot <- ggplot() + 
        geom_point(mapping = aes(x = xy, y = y), size = size) +
        geom_line(mapping = aes(x = xy, y = pred_mean), colour = "blue", linetype = "dashed") +
        geom_ribbon(mapping = aes(x = xy, 
                                  ymin = my_logistic(lb),
                                  ymax = my_logistic(ub)), fill = "blue", alpha = alpha) +
        # geom_line(mapping = aes(x = xy, y = my_logistic(lb) * m ), colour = "blue", linetype = 2) +
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        ggtitle(title)
    }
  }
  
  
  return(myplot)
  
}

## function to take data and prediction results and create a nice data frame ready for plotting
data_to_plot_df <- function(data, pred, pred_type, ci = c(0.025,0.975))
{
  plot_df <- data[data$train == FALSE,]
  # plot_df$pred <- pred$pred$pred_mean
  plot_df$estimate <- pred$pred$pred_mean
  # plot_df$pred_var <- pred$pred$pred_var
  plot_df$estimate_type <- rep("link", times = nrow(plot_df))
  
  plot_df_response <- data[data$train == FALSE,]
  plot_df_response$estimate <- pred$inverse_link(pred$pred$pred_mean)
  plot_df_response$estimate_type <- rep("response", times = nrow(plot_df_response))
  
  plot_df_lower <- data[data$train == FALSE,]
  plot_df_lower$estimate <- pred$inverse_link(pred$pred$pred_mean + qnorm(p = ci[1], mean = 0, sd = 1) * pred$pred$pred_var)
  plot_df_lower$estimate_type <- rep("lower_response", times = nrow(plot_df_response))
  
  plot_df_upper <- data[data$train == FALSE,]
  plot_df_upper$estimate <- pred$inverse_link(pred$pred$pred_mean + qnorm(p = ci[2], mean = 0, sd = 1) * pred$pred$pred_var)
  plot_df_upper$estimate_type <- rep("upper_response", times = nrow(plot_df_response))
  
  
  # plot_df$pred_response <- pred$inverse_link(pred$pred$pred_mean)
  # plot_df$lower <- pred$inverse_link(pred$pred$pred_mean + qnorm(p = ci[1], mean = 0, sd = 1) * pred$pred$pred_var)
  # plot_df$upper <- pred$inverse_link(pred$pred$pred_mean + qnorm(p = ci[2], mean = 0, sd = 1) * pred$pred$pred_var)
  plot_df_final <- rbind(plot_df, plot_df_response, plot_df_lower, plot_df_upper)
  plot_df_final$type <- rep(pred_type, times = nrow(plot_df_final))
  plot_df_final$type <- as.factor(plot_df_final$type)
  plot_df_final$estimate_type <- as.factor(plot_df_final$estimate_type)
  return(plot_df_final)
  
}
