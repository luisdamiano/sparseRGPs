## covariance functions
#' @export
cov_fun_sqrd_exp <- function(x1, x2, cov_par)
{
  sigma <- cov_par$sigma
  l <- cov_par$l

  return(sigma^2 * exp(-1/(2*l^2) * (x1 - x2)^2))
}

#' @export
mv_cov_fun_sqrd_exp <- function(x1, x2, cov_par)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l

  return(sigma^2 * exp(-1/(2*l^2) * sum((x1 - x2)^2)))
}

#' @export
cov_fun_sqrd_exp_ard <- function(x1, x2, cov_par)
{
  ## x1, x2 are row vectors
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

  return(sigma^2 * exp(-1/(2) * sum((x1 - x2)^2 / l^2)))
}

#' @export
mv_cov_fun_exp <- function(x1, x2, cov_par)
{
  ## x1, x2 are row vectors
  sigma <- cov_par$sigma
  l <- cov_par$l

  return(sigma^2 * exp(-1/(l) * sqrt(sum((x1 - x2)^2))))
}


## function to make a covariance matrix from a covariance function and input values
#' @export
make_cov_mat <- function(x, x_pred, cov_fun, cov_par, tau = 0)
{
  xfull <- rbind(x,x_pred)
  temp <- matrix(nrow = nrow(xfull), ncol = nrow(xfull))
  for(i in 1:nrow(xfull))
  {
    for(j in 1:nrow(xfull))
    {
      temp[i,j] <- cov_fun(xfull[i,], xfull[j,], cov_par) + 1*(i==j)*tau^2
    }
  }
  return(temp)
}


## matrix population function
#Rcpp::cppFunction('double cov_fun_sqrd_expC(NumericVector x1, NumericVector x2, List cov_par){
#
#double sigma = as<double>(cov_par["sigma"]);
#double l = as<double>(cov_par["l"]);
#
#return pow(sigma,2) * exp(-1/(2 * pow(l,2)) * sum(pow((x1 - x2),2)));
#
#}')

# x1 <- rnorm(n = 5)
# x2 <- rnorm(n = 5)
# sigma <- 1
# l <- 1
# cov_par_test <- list("sigma" = 1, "l" = 1)
#
# cov_fun_sqrd_expC(x1 = x1, x2 = numeric(), cov_par = cov_par_test)


## squared exponential covariance function in Rcpp IT WORKS
# Rcpp::cppFunction('NumericMatrix make_cov_matC(NumericVector x1, NumericVector x2, List cov_par, String cov_fun){
#
# if(cov_fun == "sqexp"){
#
#  if(x2.size() == 0){
#      nrow = x1.size();
#      ncol = x1.size();
#
#      NumericMatrix mat(nrow, ncol);
#
#      for(int i = 0; i < nrow; i++){
#         for(int j = 0; j < ncol; j++){
#           mat[i,j] = cov_fun_sqrd_expC(x1 = x1[i], x2 = x2[j], cov_par = cov_par);
#         }
#      }
#
#  }
#  return mat;
# }
#
# else{return NumericMatrix mat(0,0);}
#
# }')
#
# make_cov_matC(x1 = 0, x2 = 0, cov_par = cov_par_test, cov_fun = "crap")


#Rcpp::cppFunction( 'double na_fun(NumericMatrix x){
#    return x.size();
#  }'
#)
#na_fun(matrix())
