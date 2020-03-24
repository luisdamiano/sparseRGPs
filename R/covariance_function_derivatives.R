## derivatives of covariance functions

## squared exponential covariance function

## derivative wrt sigma
#' @export
dsqexp_dsigma <- function(x1, x2, cov_par, transform = FALSE)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l

  if(transform == TRUE)
  {
    dsigma_dsigmat <- sigma
    trans_fun <- function(x)
    {
      return(exp(x))
    }
    inv_trans_fun <- function(x)
    {
      return(log(x))
    }
    return(
      list("derivative" = 2 * sigma * exp( -(1/(2*l^2)) * sum((x1 - x2)^2)) * dsigma_dsigmat,
           "trans_fun" = trans_fun,
           "trans_par" = inv_trans_fun(sigma),
           "inverse_trans_fun" = inv_trans_fun)
    )
  }
  if(transform == FALSE)
  {
    return(
      list("derivative" = 2 * sigma * exp( -(1/(2*l^2)) * sum((x1 - x2)^2)))
    )
  }

}

#' @export
dsqexp_dsigma_ard <- function(x1, x2, cov_par, transform = FALSE)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma


  l <- numeric()
  for(i in 1:length(x1))
  {
    l[i] <- cov_par[[which(names(cov_par) == paste("l", i, sep = ""))]]
  }

  if(transform == TRUE)
  {
    dsigma_dsigmat <- sigma
    trans_fun <- function(x)
    {
      return(exp(x))
    }
    inv_trans_fun <- function(x)
    {
      return(log(x))
    }
    return(
      list("derivative" = 2 * sigma * exp( -(1/(2)) * sum((x1 - x2)^2 / (l^2))) * dsigma_dsigmat,
           "trans_fun" = trans_fun,
           "trans_par" = inv_trans_fun(sigma),
           "inverse_trans_fun" = inv_trans_fun)
    )
  }
  if(transform == FALSE)
  {
    return(
      list("derivative" = 2 * sigma * exp( -(1/(2*l^2)) * sum((x1 - x2)^2)))
    )
  }

}

## derivative wrt tau, the nugget sd
#' @export
dsqexp_dtau <- function(x1, x2, cov_par, transform = FALSE)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l
  tau <- cov_par$tau

  if(transform == TRUE)
  {
    dtau_dtaut <- tau
    trans_fun <- function(x)
    {
      return(exp(x))
    }
    inv_trans_fun <- function(x)
    {
      return(log(x))
    }
    return(
      list("derivative" = 2 * tau * dtau_dtaut * 1 * (all(x1 == x2)),
           "trans_fun" = trans_fun,
           "trans_par" = inv_trans_fun(tau),
           "inverse_trans_fun" = inv_trans_fun)
    )
  }
  if(transform == FALSE)
  {
    return(
      list("derivative" = 2 * tau * 1 * (all(x1 == x2)))
    )
  }

}

# dsqexp_dsigma(x1 = 1, x2 = 1, cov_par = list("sigma" = 2, "l" = 0.1), transform = FALSE)
# dsqexp_dsigma(x1 = 1, x2 = 2, cov_par = list("sigma" = 2, "l" = 0.1), transform = FALSE)


## derivative wrt l
#' @export
dsqexp_dl <- function(x1, x2, cov_par, transform = FALSE)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l

  if(transform == TRUE)
  {
    dl_dlt <- l
    trans_fun <- function(x)
    {
      return(exp(x))
    }
    inv_trans_fun <- function(x)
    {
      return(log(x))
    }
    return(
      list("derivative" = (sigma^2 * exp( (-1/(2*l^2)) * sum((x1 - x2)^2) )) *
             ( (1/(l^3)) * sum((x1 - x2)^2)) * dl_dlt,
           "trans_fun" = trans_fun,
           "trans_par" = inv_trans_fun(l),
           "inverse_trans_fun" = inv_trans_fun)
    )
  }
  if(transform == FALSE)
  {
    return(
      list("derivative" = (sigma^2 * exp( (-1/(2*l^2)) * sum((x1 - x2)^2) )) * ( (1/(l^3)) * sum((x1 - x2)^2)))
    )
  }

}

# ls <- seq(from = -1, to = 1, by = 1e-4)
# ls_plus <- ls + 1e-4 / 2
# ls_minus <- ls - 1e-4/2
# dls1 <-  numeric()
# dls2 <- numeric()
# x1 <- 1
# x2 <- 2
# sigma <- 1
# for(i in 1:length(ls))
# {
#   dls1[i] <- dsqexp_dl(x1 = x1, x2 = x2, cov_par = list("sigma" = 1, "l" = exp(ls[i])), transform = TRUE)$derivative
#   dls2[i] <- ((sigma^2 * exp( (-1/(2*exp(ls_plus[i])^2)) * sum((x1 - x2)^2) )) -
#     (sigma^2 * exp( (-1/(2*exp(ls_minus[i])^2)) * sum((x1 - x2)^2) ))) / 1e-4
# }
#
# plot(dls2, y = dls1)
# dsqexp_dl(x1 = 1, x2 = 1, cov_par = list("sigma" = 2, "l" = 0.1), transform = FALSE)
# dsqexp_dl(x1 = 1, x2 = 2, cov_par = list("sigma" = 2, "l" = 1), transform = TRUE)


## derivative wrt x2
#' @export
dsqexp_dx2 <- function(x1, x2, cov_par, transform = FALSE, bounds)
{
  ## x1, x2 are row vectors
  ## bounds is a matrix of dimension length(x2) x 2
  ##    where the left value is a lower bound and right value is an upper bound
  ##    for each dimension
  sigma <- cov_par$sigma
  l <- cov_par$l
  dx2_dx2t <- (bounds[,2] - bounds[,1]) / (( (x2 - bounds[,1] ) * (bounds[,2] - x2 )) + 1e-4)

  if(transform == TRUE)
  {
    ## function goes from real valued to bounded
    trans_fun <- function(x, bounds)
    {
      return(
          # bounds[,2] * exp(-log( (1/exp(x) + 1))) + bounds[,1] * exp(log(1 / (1 + exp(x)) + 1e-6))
        bounds[,2] * (1 / (1 + exp(-x))) + bounds[,1] * (1 / (1 + exp(x)))
        )
    }
    ## function goes from bounded to real valued
    inv_trans_fun <- function(x,bounds)
    {
      return(log((x - bounds[,1]) + 1e-4) - log((bounds[,2] - x) + 1e-4))
    }
    return(
      ifelse(test = is.matrix(x1) | is.matrix(x2),
             yes = return(
               list("derivative" = as.numeric(t((1/(l^2)) * (x1 - x2) * sigma^2 * (exp( -sum((x1 - x2)^2) / (2*l^2))))) * dx2_dx2t,
                    "trans_par" = inv_trans_fun(x = x2, bounds = bounds),
                    "trans_fun" = trans_fun,
                    "inverse_trans_fun" = inv_trans_fun)
             ),
             no = return(
               list("derivative" =  (1/(l^2)) * (x1 - x2) * sigma^2 * (exp( -sum((x1 - x2)^2) / (2*l^2))) * dx2_dx2t,
                    "trans_par" = inv_trans_fun(x = x2, bounds = bounds),
                    "trans_fun" = trans_fun,
                    "inverse_trans_fun" = inv_trans_fun)
             )
      )
    )
  }
  if(transform == FALSE)
  {
    return(
      ifelse(test = is.matrix(x1) | is.matrix(x2),
             yes = return(
               list("derivative" = as.numeric(t((1/(l^2)) * (x1 - x2) * sigma^2 * (exp( -sum((x1 - x2)^2) / (2*l^2))))),
                    "trans_par" = x2)
             ),
             no = return(
               list("derivative" =  (1/(l^2)) * (x1 - x2) * sigma^2 * (exp( -sum((x1 - x2)^2) / (2*l^2))),
                    "trans_par" = x2)
             )
      )
    )
  }
}

## derivative wrt x2
#' @export
dsqexp_dx2_ard <- function(x1, x2, cov_par, transform = FALSE, bounds)
{
  ## x1, x2 are row vectors
  ## bounds is a matrix of dimension length(x2) x 2
  ##    where the left value is a lower bound and right value is an upper bound
  ##    for each dimension
  sigma <- cov_par$sigma
  l <- numeric(length(x1))
  for(i in 1:length(x1))
  {
    l[i] <- eval(
      parse(text = eval(substitute(expr =
                                     paste("cov_par$l", b, sep = ""),
                                   env = list("b" = i)))
      )
    )
  }
  dx2_dx2t <- (bounds[,2] - bounds[,1]) / (( (x2 - bounds[,1] ) * (bounds[,2] - x2 )) + 1e-4)

  if(transform == TRUE)
  {
    ## function goes from real valued to bounded
    trans_fun <- function(x, bounds)
    {
      return(
        # bounds[,2] * exp(-log( (1/exp(x) + 1))) + bounds[,1] * exp(log(1 / (1 + exp(x)) + 1e-6))
        bounds[,2] * (1 / (1 + exp(-x))) + bounds[,1] * (1 / (1 + exp(x)))
      )
    }
    ## function goes from bounded to real valued
    inv_trans_fun <- function(x,bounds)
    {
      return(log((x - bounds[,1]) + 1e-4) - log((bounds[,2] - x) + 1e-4))
    }
    return(
      ifelse(test = is.matrix(x1) | is.matrix(x2),
             yes = return(
               list("derivative" = as.numeric(t((1/(l^2)) * (x1 - x2) * sigma^2 *
                                                  ( exp(-1/(2) * sum((x1 - x2)^2 / l^2)) ))) * dx2_dx2t,
                    "trans_par" = inv_trans_fun(x = x2, bounds = bounds),
                    "trans_fun" = trans_fun,
                    "inverse_trans_fun" = inv_trans_fun)
             ),
             no = return(
               list("derivative" =  (1/(l^2)) * (x1 - x2) * sigma^2 * (exp(-1/(2) * sum((x1 - x2)^2 / l^2))) * dx2_dx2t,
                    "trans_par" = inv_trans_fun(x = x2, bounds = bounds),
                    "trans_fun" = trans_fun,
                    "inverse_trans_fun" = inv_trans_fun)
             )
      )
    )
  }
  if(transform == FALSE)
  {
    return(
      ifelse(test = is.matrix(x1) | is.matrix(x2),
             yes = return(
               list("derivative" = as.numeric(t((1/(l^2)) * (x1 - x2) * sigma^2 * (exp(-1/(2) * sum((x1 - x2)^2 / l^2))))),
                    "trans_par" = x2)
             ),
             no = return(
               list("derivative" =  (1/(l^2)) * (x1 - x2) * sigma^2 * (exp(-1/(2) * sum((x1 - x2)^2 / l^2))),
                    "trans_par" = x2)
             )
      )
    )
  }
}

# xu <- matrix(rep(1:4, each = 1), nrow = 2, ncol = 2)
# bounds <- matrix(c(1 - 1e-3,2 + 1e-3,3 - 1e-3,4 + 1e-3), nrow = 2, ncol = 2, byrow = TRUE)
#
# dsqexp_dx2(x1 = matrix(xu[1,],nrow = 1), x2 = matrix(xu[2,], nrow = 1), cov_par = list("sigma" = 2, "l" = 2), transform = TRUE, bounds = bounds)
#
# dsqexp_dx2(x1 = matrix(c(1.25,3.25), nrow = 1),
#            x2 = matrix(c(1.75,3.75),nrow = 1),
#            cov_par = list("sigma" = 2, "l" = 2), bounds = bounds, transform = TRUE)

# blork <- mv_cov_fun_sqrd_exp(x1 = 0, x2 = -1 + 1e-4, cov_par = list("sigma" = 2, "l" = 2))
# blark <- mv_cov_fun_sqrd_exp(x1 = 0, x2 = -1 - 1e-4, cov_par = list("sigma" = 2, "l" = 2))


## EXPONENTIAL COVARIANCE FUNCTION
## squared exponential covariance function

## derivative wrt sigma
#' @export
dexp_dsigma <- function(x1, x2, cov_par, transform = FALSE)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l

  if(transform == TRUE)
  {
    dsigma_dsigmat <- sigma
    trans_fun <- function(x)
    {
      return(exp(x))
    }
    inv_trans_fun <- function(x)
    {
      return(log(x))
    }
    return(
      list("derivative" = 2 * sigma * exp( -(1/(l)) * sqrt(sum((x1 - x2)^2))) * dsigma_dsigmat,
           "trans_fun" = trans_fun,
           "trans_par" = inv_trans_fun(sigma),
           "inverse_trans_fun" = inv_trans_fun)
    )
  }
  if(transform == FALSE)
  {
    return(
      list("derivative" = 2 * sigma * exp( -(1/(l)) * sqrt(sum((x1 - x2)^2)) ) )
    )
  }

}

## derivative wrt tau, the nugget sd
#' @export
dexp_dtau <- function(x1, x2, cov_par, transform = FALSE)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l
  tau <- cov_par$tau

  if(transform == TRUE)
  {
    dtau_dtaut <- tau
    trans_fun <- function(x)
    {
      return(exp(x))
    }
    inv_trans_fun <- function(x)
    {
      return(log(x))
    }
    return(
      list("derivative" = 2 * tau * dtau_dtaut * 1 * (all(x1 == x2)),
           "trans_fun" = trans_fun,
           "trans_par" = inv_trans_fun(tau),
           "inverse_trans_fun" = inv_trans_fun)
    )
  }
  if(transform == FALSE)
  {
    return(
      list("derivative" = 2 * tau * 1 * (all(x1 == x2)))
    )
  }

}
# dsqexp_dtau(x1 = 1, x2 = 2, cov_par = list("sigma" = 2, "l" = 0.1, "tau" = 0.5), transform = TRUE)
# dsqexp_dsigma(x1 = 1, x2 = 1, cov_par = list("sigma" = 2, "l" = 0.1), transform = FALSE)
# dsqexp_dsigma(x1 = 1, x2 = 2, cov_par = list("sigma" = 2, "l" = 0.1), transform = FALSE)


## derivative wrt l
#' @export
dexp_dl <- function(x1, x2, cov_par, transform = FALSE)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l

  if(transform == TRUE)
  {
    dl_dlt <- l
    trans_fun <- function(x)
    {
      return(exp(x))
    }
    inv_trans_fun <- function(x)
    {
      return(log(x))
    }
    return(
      list("derivative" = (sigma^2 * exp( (-1/(l)) * sqrt(sum((x1 - x2)^2)) )) * ( (1/(l^2)) * sqrt(sum((x1 - x2)^2)) ) * dl_dlt,
           "trans_fun" = trans_fun,
           "trans_par" = inv_trans_fun(l),
           "inverse_trans_fun" = inv_trans_fun)
    )
  }
  if(transform == FALSE)
  {
    return(
      list("derivative" = (sigma^2 * exp( (-1/(l)) * sqrt(sum((x1 - x2)^2)) )) * ( (1/(l^2)) * sqrt(sum((x1 - x2)^2))) )
    )
  }

}


