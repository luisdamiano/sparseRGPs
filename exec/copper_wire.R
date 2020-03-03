## Code to replicate the simulated copper wire results

## uncomment this next line of code if you used
##    devtools::install_github("nategarton13/sparseRGPs") to download the
##    R package
#library(sparseRGPs)

##  comment this next line of code if you used
##    devtools::install_github("nategarton13/sparseRGPs") to download the
##    R package
devtools::load_all()
library(ggplot2)
library(reshape2)
library(tidyr)
library(gridExtra)
library(caret)
library(dplyr)
library(mvtnorm)

## dissimilarity function
delta_fun <- function(x,y)
{
  return(sum((x - y)^2))
}

## measures dissimilarity between e_k and e_u
##    eagg <- rbind(e_u, e_k, e_3, ..., e_I)
delta_wrapper <- function(eagg)
{
  eagg <- t(eagg)
  eu <- eagg[1,]
  ek <- eagg[2,]

  ek_delta <- delta_fun(x = ek, y = eu)
  return(ek_delta)
}

## score functions
##    eagg <- rbind(e_u, e_k, e_3, ..., e_I)
smin_fun <- function(eagg)
{
  # eu is the vector of measurements on the unknown source piece of evidence
  # ek is the vector of measurements on the known source piece of evidence
  #   that we care about matching with eu
  # e is a matrix with rows equal to the sources
  #   columns are variables with order matching eu
  eagg <- t(eagg)
  eu <- eagg[1,]
  ek <- eagg[2,]
  e <- eagg[3:nrow(eagg),]

  e_delta <- apply(X = e, MARGIN = 1, FUN = delta_fun, y = eu)
  ek_delta <- delta_fun(x = ek, y = eu)

  return(log(ek_delta) - log(min(e_delta)))
}

##    eagg <- rbind(e_u, e_k, e_3, ..., e_I)
smax_fun <- function(eagg)
{
  # eu is the vector of measurements on the unknown source piece of evidence
  # ek is the vector of measurements on the known source piece of evidence
  #   that we care about matching with eu
  # e is a matrix with rows equal to the sources
  #   columns are variables with order matching eu
  eagg <- t(eagg)
  eu <- eagg[1,]
  ek <- eagg[2,]
  e <- eagg[3:nrow(eagg),]

  e_delta <- apply(X = e, MARGIN = 1, FUN = delta_fun, y = eu)
  ek_delta <- delta_fun(x = ek, y = eu)

  return(log(ek_delta) - log(max(e_delta)))
}

##    eagg <- rbind(e_u, e_k, e_3, ..., e_I)
savg_fun <- function(eagg)
{
  # eu is the vector of measurements on the unknown source piece of evidence
  # ek is the vector of measurements on the known source piece of evidence
  #   that we care about matching with eu
  # e is a matrix with rows equal to the sources
  #   columns are variables with order matching eu
  eagg <- t(eagg)
  eu <- eagg[1,]
  ek <- eagg[2,]
  e <- eagg[3:nrow(eagg),]

  e_delta <- apply(X = e, MARGIN = 1, FUN = delta_fun, y = eu)
  ek_delta <- delta_fun(x = ek, y = eu)

  return(log(ek_delta) - log(mean(e_delta)))
}

## group means
## (for Ag, Sb, Pb, Bi, Co, Ni, As, and Se, respectively)
mu1 <- c(7.429, 0.761, 0.774, 0.070, 0.030, 2.033, 0.917, 0.756)
mu2 <- c(4.524, 0.465, 0.687, 0.037, 0.28, 1.065, 0.512, 0.536)

## within source covariance matrix
data("within_covariance")

## between source covariance matrix
data("between_covariance")

## probability of class 1 membership
pi1 <- 0.172

## generate source function
##    This function generates mean vectors for a given number of sources
##    n: the number of sources
##    pi1: the probability of group one in the between source Gaussian
##        mixture model
##    mu1: the mean of group one in the between source Gaussian mixture model
##    mu2: the mean of group two in the between source Gaussian mixture model
##    between_cov: the between source covariance matrix in the between source
##        Gaussian mixture model
generate_source <- function(n, pi1, mu1, mu2, between_cov)
{
  # sample the group
  group <- rbinom(n = n, size = 1, prob = pi1)
  group[group == 0] <- 2

  means <- matrix(nrow = n, ncol = length(mu1))
  samples <- matrix(nrow = n, ncol = length(mu1))
  for(i in 1:n)
  {
    ## generate the mean for the right group
    means[i,] <- 1*(group[i] == 1)*mu1 + 1*(group[i] == 2)*mu2

    ## sample from the correct gaussian to get the source mean
    samples[i,] <- mvtnorm::rmvnorm(n = 1, mean = means[i,], sigma = between_cov)
  }

  return(samples)

}

## feature-based likelihood ratio function
##  eu: vector of the unknown source chemical concentrations
##  mu_mat: rbind(mu_k, mu_2, ..., e_{I - 1}), the matrix of
##    mean chemical concentrations of the known sources
##    within_cov: the within source covariance matrix
##  probs: the mixing probabilities for H_{ss} = d
##  log: if TRUE, return the log-likelihood ratio
lr_fun <- function(eu, mu_mat, within_cov, probs = NA, log = TRUE)
{
  if(is.na(probs))
  {
    probs <- rep(1/(nrow(mu_mat) - 1), times = nrow(mu_mat) - 1)
  }
  if(log == TRUE)
  {
    numer <- mvtnorm::dmvnorm(x = eu, mean = mu_mat[1,], sigma = within_cov, log = log)
    denom <- log(sum(apply(X = mu_mat[-1,], MARGIN = 1, FUN = function(mu, sigma, eu)
    {
      return(
        mvtnorm::dmvnorm(x = eu, mean = mu, sigma = sigma, log = FALSE)
      )
    }, eu = eu, sigma = within_cov) * probs))

    return(numer - denom)
  }
  if(log == FALSE)
  {
    numer <- mvtnorm::dmvnorm(x = eu, mean = mu_mat[1,], sigma = within_cov, log = log)
    denom <- (sum(apply(X = mu_mat[-1,], MARGIN = 1, FUN = function(mu, sigma, eu)
    {
      return(
        mvtnorm::dmvnorm(x = eu, mean = mu, sigma = sigma, log = FALSE)
      )
    }, sigma = within_cov, eu = eu) * probs))

    return(numer / denom)
  }
}

## generate data under prosecutions hypothesis
##    n: the number of samples to generate
##    mue: the matrix of mean vectors for the sources in a database, i.e.
##      not e_u or e_k
##    muek: the mean vector of chemical concentrations of the known source
##      evidence which has the same source as the unknown source evidence
##    mueu: the mean vector of chemical concentrations of the unknown source
##      evidence which has the same source as the known source evidence
##      (technically this argument is redundant)
generate_hp <- function(n, mue, muek, mueu, within_cov)
{
  ## store samples in 3D array
  samples_array <- array(dim = c(n, length(muek), nrow(mue) + 2))

  ## sample reps for each source individually
  ## first dimension 3 is unknown source evidence
  samples_array[,,1] <- (mvtnorm::rmvnorm(n = n, mean = mueu, sigma = within_cov))

  ## second dimension 3 is known source evidence
  samples_array[,,2] <- (mvtnorm::rmvnorm(n = n, mean = muek, sigma = within_cov))

  for(i in 3:(nrow(mue) + 2))
  {
    samples_array[,,i] <- (mvtnorm::rmvnorm(n = n, mean = mue[i - 2,], sigma = within_cov))
  }

  return(samples_array)

}

## generate data under defense hypothesis
##    n: the number of samples to generate
##    mue: the matrix of mean vectors of chemical concentrations of the
##      database sources, note that the mean of the unknown source evidence
##      comes from one of these sources under H_{ss} = d
##    muek: the mean vector of chemical concentrations of the known source
##      evidence which, under H_{ss} = p, has the same source as
##      the known source evidence
##    probs: the mixing probabilities for H_{ss} = d
##    within_cov: the within source covariance matrix
generate_hd <- function(n, mue, muek, within_cov, probs = NA)
{
  if(is.na(probs))
  {
    probs <- rep(1/nrow(mue), times = nrow(mue))
  }

  ## probs is a vector of prior probabilities over the sources mue
  source_eu <- sample.int(n = nrow(mue), size = n, replace = TRUE, prob = probs)

  ## store samples in 3D array
  samples_array <- array(dim = c(n, length(muek), nrow(mue) + 2))

  ## sample reps for each source individually
  ## first dimension 3 is unknown source evidence
  meanmat_eu <- mue[source_eu,]
  # samples_array[,,1] <- mvtnorm::rmvnorm(n = 1, mean = meanmat_eu, sigma = within_cov)
  samples_array[,,1] <- t(apply(X = meanmat_eu, MARGIN = 1, FUN =
                                  function(mu, sigma){mvtnorm::rmvnorm(n = 1, mean = mu, sigma = sigma)},
                                sigma = within_cov))

  ## second dimension 3 is known source evidence
  samples_array[,,2] <- (mvtnorm::rmvnorm(n = n, mean = muek, sigma = within_cov))

  for(i in 3:(nrow(mue) + 2))
  {
    samples_array[,,i] <- (mvtnorm::rmvnorm(n = n, mean = mue[i - 2,], sigma = within_cov))
  }

  return(samples_array)

}

## perform several runs
runs <- 3 ## to replicate exact paper results, set this to 5
n <- 200 ## to replicate exact paper results, set this to 2000
N_A <- 500 ## number of known sources
I <- N_A + 1 ## number of total pieces of evidence
TTmax <- 40 ## algorithmic parameters
TTmin <- 5 ## algorithmic parameters
maxit <- 500 ## algorithmic parameters
maxit_nr <- 200 ## algorithmic parameters
maxknot <- 75 ## algorithmic parameters
delta <- 1e-2 ## algorithmic parameters
epsilon <- 1e-4 ## algorithmic parameters
results_table <- data.frame("score" = "temp",
                            "hypothesis" = "temp",
                            "true-KL" = 0,
                            "true-KL-SD" = 0,
                            "score-KL" = 0,
                            "score-KL-SD" = 0,
                            "RMSE" = 0, "Run" = 0)

system.time(for(i in 1:runs)
{
  set.seed(1307 + i)
  source_prior <- rep(1/(N_A - 1), times = N_A - 1)

  ## generate sources
  sources <- generate_source(n = N_A,
                             pi1 = pi1,
                             mu1 = mu1,
                             mu2 = mu2,
                             between_cov = between_covariance)
  muek <- sources[1,]
  mue <- sources[2:N_A,]

  ## generate data
  train_data_hp <- generate_hp(n = n,
                               mue = mue,
                               muek = muek,
                               mueu = muek,
                               within_cov = within_covariance)
  train_data_hd <- generate_hd(n = n,
                               mue = mue,
                               muek = muek,
                               within_cov = within_covariance,
                               probs = source_prior)

  ## get scores
  sdelta_hp <- apply(X = train_data_hp, MARGIN = 1, FUN = delta_wrapper)
  smin_hp <- apply(X = train_data_hp, MARGIN = 1, FUN = smin_fun)
  savg_hp <- apply(X = train_data_hp, MARGIN = 1, FUN = savg_fun)
  smax_hp <- apply(X = train_data_hp, MARGIN = 1, FUN = smax_fun)

  sdelta_hd <- apply(X = train_data_hd, MARGIN = 1, FUN = delta_wrapper)
  smin_hd <- apply(X = train_data_hd, MARGIN = 1, FUN = smin_fun)
  savg_hd <- apply(X = train_data_hd, MARGIN = 1, FUN = savg_fun)
  smax_hd <- apply(X = train_data_hd, MARGIN = 1, FUN = smax_fun)

  ## generate test data
  test_data_hp <- generate_hp(n = n,
                              mue = mue,
                              muek = muek,
                              mueu = muek,
                              within_cov = within_covariance)
  test_data_hd <- generate_hd(n = n,
                              mue = mue,
                              muek = muek,
                              within_cov = within_covariance,
                              probs = source_prior)

  ## get scores
  sdelta_hp_test <- apply(X = test_data_hp, MARGIN = 1, FUN = delta_wrapper)
  smin_hp_test <- apply(X = test_data_hp, MARGIN = 1, FUN = smin_fun)
  savg_hp_test <- apply(X = test_data_hp, MARGIN = 1, FUN = savg_fun)
  smax_hp_test <- apply(X = test_data_hp, MARGIN = 1, FUN = smax_fun)

  sdelta_hd_test <- apply(X = test_data_hd, MARGIN = 1, FUN = delta_wrapper)
  smin_hd_test <- apply(X = test_data_hd, MARGIN = 1, FUN = smin_fun)
  savg_hd_test <- apply(X = test_data_hd, MARGIN = 1, FUN = savg_fun)
  smax_hd_test <- apply(X = test_data_hd, MARGIN = 1, FUN = smax_fun)

  y_test <- rep(c(1,0), each = n)


  ## actual LR
  log_lr_hp <- apply(X = test_data_hp[,,1],
                     MARGIN = 1,
                     FUN = lr_fun,
                     mu_mat = sources,
                     within_cov = within_covariance,
                     probs = NA, log = TRUE)

  log_lr_hd <- apply(X = test_data_hd[,,1],
                     MARGIN = 1,
                     FUN = lr_fun,
                     mu_mat = sources,
                     within_cov = within_covariance,
                     probs = NA, log = TRUE)

  ########################################################
  ## train sparse GPs to get SLRs
  ########################################################
  sdelta <- log(c(sdelta_hp, sdelta_hd))
  y <- rep(c(1,0), each = n)
  xu0_delta <- kmeans(x = sdelta[1:n], centers = 5)$centers
  xu1_delta <- kmeans(x = sdelta[(n + 1):(2*n)], centers = 5)$centers
  xu_delta <- rbind(xu0_delta, xu1_delta)
  cp_start <- list("sigma" = 10, "l" = 6, "tau" = 1e-3)

  ## The "delta" score
  set.seed(1308)
  system.time(gp_mod_sdelta <- optimize_gp(y = y,
                                           xy = sdelta, cov_fun = "sqexp",
                                           cov_par_start = cp_start,
                                           mu = rep(0, times = length(y)), family = "bernoulli",
                                           nugget = TRUE, sparse = TRUE, xu_opt = "oat",
                                           xu = xu_delta, muu = rep(0, times = nrow(xu_delta)),
                                           opt = list("Tmin" = TTmin,
                                                      "Tmax" = TTmax,
                                                      "delta" = delta,
                                                      "obj_tol" = 1e-3,
                                                      "maxknot" = maxknot,
                                                      "maxit" = maxit,
                                                      "maxit_nr" = maxit_nr,
                                                      epsilon = epsilon,
                                                      "grad_tol" = Inf, delta = delta),
                                           verbose = FALSE, file_path = NULL))

  print("delta mod done")

  sdelta_test <- log(c(sdelta_hp_test, sdelta_hd_test))
  pred_gp_sdelta <- predict_gp(mod = gp_mod_sdelta,
                               x_pred = sdelta_test,
                               mu_pred = rep(0, times = length(y_test)),
                               full_cov = FALSE, vi = FALSE)

  sgp_hp_test_sdelta <- pred_gp_sdelta$pred$pred_mean[1:n]
  sgp_hd_test_sdelta <- pred_gp_sdelta$pred$pred_mean[(n + 1):(2*n)]

  prior_p <- mean(y)

  log_slr_gp_sdelta <- log(my_logistic(pred_gp_sdelta$pred$pred_mean)) -
    log(1 - my_logistic(pred_gp_sdelta$pred$pred_mean)) -
    log(prior_p / (1 - prior_p))

  ## "min" score
  smin <- c(smin_hp, smin_hd)
  y <- rep(c(1,0), each = n)
  xu0_min <- kmeans(x = smin[1:n], centers = 5)$centers
  xu1_min <- kmeans(x = smin[(n + 1):(2*n)], centers = 5)$centers
  xu_min <- rbind(xu0_min, xu1_min)
  cp_start <- list("sigma" = 10, "l" = 6, "tau" = 1e-3)

  set.seed(1308)
  system.time(gp_mod_smin <- optimize_gp(y = y,
                                         xy = smin, cov_fun = "sqexp",
                                         cov_par_start = cp_start,
                                         mu = rep(0, times = length(y)), family = "bernoulli",
                                         nugget = TRUE, sparse = TRUE, xu_opt = "oat",
                                         xu = xu_min, muu = rep(0, times = nrow(xu_min)),
                                         opt = list("Tmin" = TTmin,
                                                    "Tmax" = TTmax,
                                                    "delta" = delta,
                                                    "obj_tol" = 1e-3,
                                                    "maxknot" = maxknot,
                                                    "maxit" = maxit,
                                                    "maxit_nr" = maxit_nr,
                                                    epsilon = epsilon,
                                                    "grad_tol" = Inf, delta = delta),
                                         verbose = FALSE, file_path = NULL))

  print("min mod done")


  smin_test <- c(smin_hp_test, smin_hd_test)
  pred_gp_smin <- predict_gp(mod = gp_mod_smin,
                             x_pred = smin_test,
                             mu_pred = rep(0, times = length(y_test)),
                             full_cov = FALSE, vi = FALSE)

  sgp_hp_test_smin <- pred_gp_smin$pred$pred_mean[1:n]
  sgp_hd_test_smin <- pred_gp_smin$pred$pred_mean[(n + 1):(2*n)]

  prior_p <- mean(y)

  log_slr_gp_smin <- log(my_logistic(pred_gp_smin$pred$pred_mean)) -
    log(1 - my_logistic(pred_gp_smin$pred$pred_mean)) -
    log(prior_p / (1 - prior_p))


  ## the "average" score
  savg <- c(savg_hp, savg_hd)
  y <- rep(c(1,0), each = n)
  # xu <- xy[sample.int(n = nrow(xy), size = 5, replace = FALSE),]
  xu0_avg <- kmeans(x = savg[1:n], centers = 5)$centers
  xu1_avg <- kmeans(x = savg[(n + 1):(2*n)], centers = 5)$centers
  xu_avg <- rbind(xu0_avg, xu1_avg)
  cp_start <- list("sigma" = 10, "l" = 6, "tau" = 1e-3)

  set.seed(1308)
  system.time(gp_mod_savg <- optimize_gp(y = y,
                                         xy = savg, cov_fun = "sqexp",
                                         cov_par_start = cp_start,
                                         mu = rep(0, times = length(y)), family = "bernoulli",
                                         nugget = TRUE, sparse = TRUE, xu_opt = "oat",
                                         xu = xu_avg, muu = rep(0, times = nrow(xu_avg)),
                                         opt = list("Tmin" = TTmin,
                                                    "Tmax" = TTmax,
                                                    "delta" = delta,
                                                    "obj_tol" = 1e-3,
                                                    "maxknot" = maxknot,
                                                    "maxit" = maxit,
                                                    "maxit_nr" = maxit_nr,
                                                    epsilon = epsilon,
                                                    "grad_tol" = Inf, delta = delta),
                                         verbose = FALSE, file_path = NULL))

  print("avg mod done")


  savg_test <- c(savg_hp_test, savg_hd_test)
  pred_gp_savg <- predict_gp(mod = gp_mod_savg,
                             x_pred = savg_test,
                             mu_pred = rep(0, times = length(y_test)),
                             full_cov = FALSE, vi = FALSE)

  sgp_hp_test_savg <- pred_gp_savg$pred$pred_mean[1:n]
  sgp_hd_test_savg <- pred_gp_savg$pred$pred_mean[(n + 1):(2*n)]

  prior_p <- mean(y)

  log_slr_gp_savg <- log(my_logistic(pred_gp_savg$pred$pred_mean)) -
    log(1 - my_logistic(pred_gp_savg$pred$pred_mean)) -
    log(prior_p / (1 - prior_p))

  ## the "max" score
  smax <- c(smax_hp, smax_hd)
  y <- rep(c(1,0), each = n)
  # xu <- xy[sample.int(n = nrow(xy), size = 5, replace = FALSE),]
  xu0_max <- kmeans(x = smax[1:n], centers = 5)$centers
  xu1_max <- kmeans(x = smax[(n + 1):(2*n)], centers = 5)$centers
  xu_max <- rbind(xu0_max, xu1_max)
  cp_start <- list("sigma" = 10, "l" = 6, "tau" = 1e-3)

  # maxknot <- 20

  # fpmax <- "RA2019/RA2019code/sparse_gp_smax_ex1.rds"
  set.seed(1308)
  system.time(gp_mod_smax <- optimize_gp(y = y,
                                         xy = smax, cov_fun = "sqexp",
                                         cov_par_start = cp_start,
                                         mu = rep(0, times = length(y)), family = "bernoulli",
                                         nugget = TRUE, sparse = TRUE, xu_opt = "oat",
                                         xu = xu_max, muu = rep(0, times = nrow(xu_max)),
                                         opt = list("Tmin" = TTmin,
                                                    "Tmax" = TTmax,
                                                    "delta" = delta,
                                                    "obj_tol" = 1e-3,
                                                    "maxknot" = maxknot,
                                                    "maxit" = maxit,
                                                    "maxit_nr" = maxit_nr,
                                                    epsilon = epsilon,
                                                    "grad_tol" = Inf, delta = delta),
                                         verbose = FALSE, file_path = NULL))

  print("max mod done")


  smax_test <- c(smax_hp_test, smax_hd_test)
  pred_gp_smax <- predict_gp(mod = gp_mod_smax,
                             x_pred = smax_test,
                             mu_pred = rep(0, times = length(y_test)),
                             full_cov = FALSE, vi = FALSE)

  sgp_hp_test_smax <- pred_gp_smax$pred$pred_mean[1:n]
  sgp_hd_test_smax <- pred_gp_smax$pred$pred_mean[(n + 1):(2*n)]

  log_slr_gp_smax <- log(my_logistic(pred_gp_smax$pred$pred_mean)) -
    log(1 - my_logistic(pred_gp_smax$pred$pred_mean)) -
    log(prior_p / (1 - prior_p))

  ## the "stacked" score
  xy <- cbind(c(smin_hp, smin_hd), c(savg_hp, savg_hd), c(smax_hp, smax_hd))
  L <- t(chol(x = cov(xy)))

  ## transform predictors to be uncorrelated
  xy_trans <- t(solve(a = L, b = t(xy)))

  y <- rep(c(1,0), each = n)
  # xu <- xy[sample.int(n = nrow(xy), size = 5, replace = FALSE),]
  xu0 <- kmeans(x = xy_trans[1:n,], centers = 5)$centers
  xu1 <- kmeans(x = xy_trans[(n + 1):(2*n),], centers = 5)$centers
  xu <- rbind(xu0, xu1)
  cp_start_ard <- list("sigma" = 10, "l1" = 2, "l2" = 2, "l3" = 2, "tau" = 1e-3)

  system.time(gp_mod <- optimize_gp(y = y,
                                    xy = xy_trans, cov_fun = "ard",
                                    cov_par_start = cp_start_ard,
                                    mu = rep(0, times = nrow(xy)), family = "bernoulli",
                                    nugget = TRUE, sparse = TRUE, xu_opt = "oat",
                                    xu = xu, muu = rep(0, times = nrow(xu)),
                                    opt = list("Tmin" = TTmin,
                                               "Tmax" = TTmax,
                                               "delta" = delta,
                                               "obj_tol" = 1e-3,
                                               "maxknot" = maxknot,
                                               "maxit" = maxit,
                                               maxit_nr = maxit_nr,
                                               "epsilon" = epsilon,
                                               "grad_tol" = Inf, delta = delta),
                                    verbose = FALSE, file_path = NULL))

  print(dim(gp_mod$results$xu))
  print("agg mod done")

  ## get gp predictions
  xy_test <- cbind(c(smin_hp_test, smin_hd_test), c(savg_hp_test, savg_hd_test),
                   c(smax_hp_test, smax_hd_test))

  xy_test_trans <- t(solve(a = L, b = t(xy_test)))
  y_test <- rep(c(1,0), each = n)

  pred_gp <- predict_gp(mod = gp_mod,
                        x_pred = xy_test_trans,
                        mu_pred = rep(0, times = nrow(xy_test)),
                        full_cov = FALSE, vi = FALSE)

  sgp_hp_test <- pred_gp$pred$pred_mean[1:n]
  sgp_hd_test <- pred_gp$pred$pred_mean[(n + 1):(2*n)]

  prior_p <- 1/2

  log_slr_gp <- log(my_logistic(pred_gp$pred$pred_mean)) -
    log(1 - my_logistic(pred_gp$pred$pred_mean)) -
    log(prior_p / (1 - prior_p))

  ## slr data frame
  log_lr <- c(log_lr_hp, log_lr_hd)
  lr_df <- data.frame("lr" = c(log_lr, log_slr_gp_sdelta, log_slr_gp_smin, log_slr_gp_savg,
                               log_slr_gp_smax, log_slr_gp),
                      "type" = c(rep(c("actual", "delta", "min", "avg", "max", "gp"), each = 2*n) ),
                      "hypothesis" = rep(rep(c("P","D"), each = n), times = 6))

  rmse_slr_delta_p <- sqrt(mean((lr_df[lr_df$hypothesis == "P" & lr_df$type == "delta",]$lr -
                                   lr_df[lr_df$hypothesis == "P" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_delta_d <- sqrt(mean((lr_df[lr_df$hypothesis == "D" & lr_df$type == "delta",]$lr -
                                   lr_df[lr_df$hypothesis == "D" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_min_p <- sqrt(mean((lr_df[lr_df$hypothesis == "P" & lr_df$type == "min",]$lr -
                                 lr_df[lr_df$hypothesis == "P" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_min_d <- sqrt(mean((lr_df[lr_df$hypothesis == "D" & lr_df$type == "min",]$lr -
                                 lr_df[lr_df$hypothesis == "D" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_avg_p <- sqrt(mean((lr_df[lr_df$hypothesis == "P" & lr_df$type == "avg",]$lr -
                                 lr_df[lr_df$hypothesis == "P" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_avg_d <- sqrt(mean((lr_df[lr_df$hypothesis == "D" & lr_df$type == "avg",]$lr -
                                 lr_df[lr_df$hypothesis == "D" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_max_p <- sqrt(mean((lr_df[lr_df$hypothesis == "P" & lr_df$type == "max",]$lr -
                                 lr_df[lr_df$hypothesis == "P" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_max_d <- sqrt(mean((lr_df[lr_df$hypothesis == "D" & lr_df$type == "max",]$lr -
                                 lr_df[lr_df$hypothesis == "D" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_gp_p <- sqrt(mean((lr_df[lr_df$hypothesis == "P" & lr_df$type == "gp",]$lr -
                                lr_df[lr_df$hypothesis == "P" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))
  rmse_slr_gp_d <- sqrt(mean((lr_df[lr_df$hypothesis == "D" & lr_df$type == "gp",]$lr -
                                lr_df[lr_df$hypothesis == "D" & lr_df$type == "actual",]$lr)^2, na.rm = TRUE))



  kl_p <- mean(lr_df[lr_df$hypothesis == "P" & lr_df$type == "actual",]$lr)
  kl_d <- mean(-lr_df[lr_df$hypothesis == "D" & lr_df$type == "actual",]$lr)
  sd_kl_p <- sd(lr_df[lr_df$hypothesis == "P" & lr_df$type == "actual",]$lr) / sqrt(n)
  sd_kl_d <- sd(lr_df[lr_df$hypothesis == "D" & lr_df$type == "actual",]$lr) / sqrt(n)

  kl_min_p <- mean(lr_df[lr_df$hypothesis == "P" & lr_df$type == "min",]$lr)
  kl_min_d <- mean(-lr_df[lr_df$hypothesis == "D" & lr_df$type == "min",]$lr)
  sd_kl_min_p <- sd(lr_df[lr_df$hypothesis == "P" & lr_df$type == "min",]$lr) / sqrt(n)
  sd_kl_min_d <- sd(lr_df[lr_df$hypothesis == "D" & lr_df$type == "min",]$lr) / sqrt(n)

  kl_avg_p <- mean(lr_df[lr_df$hypothesis == "P" & lr_df$type == "avg",]$lr)
  kl_avg_d <- mean(-lr_df[lr_df$hypothesis == "D" & lr_df$type == "avg",]$lr)
  sd_kl_avg_p <- sd(lr_df[lr_df$hypothesis == "P" & lr_df$type == "avg",]$lr) / sqrt(n)
  sd_kl_avg_d <- sd(lr_df[lr_df$hypothesis == "D" & lr_df$type == "avg",]$lr) / sqrt(n)

  kl_max_p <- mean(lr_df[lr_df$hypothesis == "P" & lr_df$type == "max",]$lr)
  kl_max_d <- mean(-lr_df[lr_df$hypothesis == "D" & lr_df$type == "max",]$lr)
  sd_kl_max_p <- sd(lr_df[lr_df$hypothesis == "P" & lr_df$type == "max",]$lr) / sqrt(n)
  sd_kl_max_d <- sd(lr_df[lr_df$hypothesis == "D" & lr_df$type == "max",]$lr) / sqrt(n)

  kl_gp_p <- mean(lr_df[lr_df$hypothesis == "P" & lr_df$type == "gp",]$lr)
  kl_gp_d <- mean(-lr_df[lr_df$hypothesis == "D" & lr_df$type == "gp",]$lr)
  sd_kl_gp_p <- sd(lr_df[lr_df$hypothesis == "P" & lr_df$type == "gp",]$lr) / sqrt(n)
  sd_kl_gp_d <- sd(lr_df[lr_df$hypothesis == "D" & lr_df$type == "gp",]$lr) / sqrt(n)

  kl_delta_p <- mean(lr_df[lr_df$hypothesis == "P" & lr_df$type == "delta",]$lr)
  kl_delta_d <- mean(-lr_df[lr_df$hypothesis == "D" & lr_df$type == "delta",]$lr)
  sd_kl_delta_p <- sd(lr_df[lr_df$hypothesis == "P" & lr_df$type == "delta",]$lr) / sqrt(n)
  sd_kl_delta_d <- sd(lr_df[lr_df$hypothesis == "D" & lr_df$type == "delta",]$lr) / sqrt(n)

  temp_df <- data.frame("score" = rep(c("delta", "min", "avg", "max", "gp"), times = 2),
                        "hypothesis" = rep(c("P","D"), each = 5),
                        "true-KL" = rep(c(kl_p, kl_d), each = 5),
                        "true-KL-SD" = rep(c(sd_kl_p, sd_kl_d), each = 5),
                        "score-KL" = c(kl_delta_p, kl_min_p, kl_avg_p, kl_max_p, kl_gp_p,
                                       kl_delta_d, kl_min_d, kl_avg_d, kl_max_d, kl_gp_d),
                        "score-KL-SD" = c(sd_kl_delta_p, sd_kl_min_p, sd_kl_avg_p, sd_kl_max_p, sd_kl_gp_p,
                                          sd_kl_delta_d, sd_kl_min_d, sd_kl_avg_d, sd_kl_max_d, sd_kl_gp_d),
                        "RMSE" = c(rmse_slr_delta_p, rmse_slr_min_p, rmse_slr_avg_p, rmse_slr_max_p, rmse_slr_gp_p,
                                   rmse_slr_delta_d, rmse_slr_min_d, rmse_slr_avg_d, rmse_slr_max_d, rmse_slr_gp_d),
                        "Run" = rep(i, times = 2*5))

  results_table <- rbind(results_table, temp_df)

  print(paste("run " , i , " complete", sep = " "))

})

# results_table
# write.csv(x = results_table[-1,], file = "/home/nmgarton/SLR/RA2019/RA2019code/copper_wire_results_NA500.csv")

results_table
