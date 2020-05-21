#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cov_fun_sqrd_expC(NumericVector x1, NumericVector x2, List cov_par){

double sigma = as<double>(cov_par["sigma"]);
double l = as<double>(cov_par["l"]);

return pow(sigma,2) * exp(-1/(2 * pow(l,2)) * sum(pow((x1 - x2),2)));

}


// [[Rcpp::export]]
double cov_fun_sqrd_exp_ardC(NumericVector x1, NumericVector x2,
                             List cov_par,
                             StringVector lnames){

  double sigma = as<double>(cov_par["sigma"]);
  int d = x1.size();
  NumericVector l(d);

  for(int i = 0; i < d; i++)
  {
    String temp = lnames[i];
    l[i] = as<double>(cov_par[temp]);
  }

  //Rcout << l << "\n";
  //NumericVector jim(d);
  //jim = -pow((x1 - x2) / l , 2) / 2;
  //double bob;
  //bob = exp(-1/2 * sum(pow((x1 - x2) / l , 2)));
  //Rcout << jim << "\n";

  //Rf_PrintValue(pow((x1 - x2) / l ,2));

  //l = as<double>(cov_par["l"]);
  return pow(sigma,2) * exp(-sum(pow((x1 - x2) / l , 2)) / 2);

}

// [[Rcpp::export]]
double cov_fun_expC(NumericVector x1, NumericVector x2, List cov_par){

  double sigma = as<double>(cov_par["sigma"]);
  double l = as<double>(cov_par["l"]);

  return pow(sigma,2) * exp(-1/l * sum(abs(x1 - x2)));

}

//' Make a covariance matrix with  covariance function
//'
//' @param x Matrix of observed data values
//' @param x_pred Potentially an empty matrix. If not empty, this returns
//' the covariance between x and x_pred.
//' @param cov_par Named list of covariance parameters. Need
//' "sigma", "tau", and "l1", ..., "ld" where d is the input dimension
//' @param cov_fun A character string specifying the covariance function.
//' Currently only works with "sqexp".
//' @param delta Scalar value that functions as a fixed nugget used to
//' stabilize matrix inverses and decompositions.
//' @return Covariance matrix at the observed data locations if
//' x_pred = matrix(), or the covariances between x and x_pred if
//' x_pred is not an empty matrix.
//' @useDynLib sparseRGPs
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericMatrix make_cov_matC(NumericMatrix x,
                            NumericMatrix x_pred,
                            List cov_par,
                            String cov_fun,
                            double delta){
  // select the squared exponential covariance function
  if(cov_fun == "sqexp"){

    // if the length of x_pred is zero, make variance covariance matrix for x
    if(NumericMatrix::is_na(x_pred(0,0)) == 1){
      int nrow = x.nrow();
      int ncol = x.nrow();

      NumericMatrix mat(nrow, ncol);

      for(int i = 0; i < nrow; i++){

        for(int j = 0; j < ncol; j++){
          if(i == j){
            mat(i , j) = cov_fun_sqrd_expC(x(i, _), x(j, _), cov_par) + pow(as<double>(cov_par["tau"]),2) + delta;
          }
          else{
            mat(i , j) = cov_fun_sqrd_expC(x(i, _), x(j, _), cov_par);
          }
        }

      }
      return mat;
    }
    // if the length of x_pred is nonzero, make the covariance matrix where x is rows, x_pred is columns
    if(NumericMatrix::is_na(x_pred(0,0)) != 1){
      int nrow = x.nrow();
      int ncol = x_pred.nrow();

      NumericMatrix mat(nrow, ncol);

      for(int i = 0; i < nrow; i++){

        for(int j = 0; j < ncol; j++){
          mat(i , j) = cov_fun_sqrd_expC(x(i, _),x_pred(j,_), cov_par);
        }

      }
      return mat;
    }
  }

  // select the exponential covariance function
  if(cov_fun == "exp"){

    // if the length of x_pred is zero, make variance covariance matrix for x
    if(NumericMatrix::is_na(x_pred(0,0)) == 1){
      int nrow = x.nrow();
      int ncol = x.nrow();

      NumericMatrix mat(nrow, ncol);

      for(int i = 0; i < nrow; i++){

        for(int j = 0; j < ncol; j++){
          if(i == j){
            mat(i , j) = cov_fun_expC(x(i, _), x(j, _), cov_par) + pow(as<double>(cov_par["tau"]),2) + delta;
          }
          else{
            mat(i , j) = cov_fun_expC(x(i, _), x(j, _), cov_par);
          }
        }

      }
      return mat;
    }
    // if the length of x_pred is nonzero, make the covariance matrix where x is rows, x_pred is columns
    if(NumericMatrix::is_na(x_pred(0,0)) != 1){
      int nrow = x.nrow();
      int ncol = x_pred.nrow();

      NumericMatrix mat(nrow, ncol);

      for(int i = 0; i < nrow; i++){

        for(int j = 0; j < ncol; j++){
          mat(i , j) = cov_fun_expC(x(i, _),x_pred(j,_), cov_par);
        }

      }
      return mat;
    }
  }

  else{
    NumericMatrix mat(0,0);
    Rcerr << "Error: invalid covariance function";
    return mat;
    }
  NumericMatrix mat(0,0);
  Rcerr << "Error: invalid covariance function";
  return mat;
}

//' Make a covariance matrix with  covariance function
//'
//' @param x Matrix of observed data values
//' @param x_pred Potentially an empty matrix. If not empty, this returns
//' the covariance between x and x_pred.
//' @param cov_par Named list of covariance parameters. Need
//' "sigma", "tau", and "l1", ..., "ld" where d is the input dimension
//' @param cov_fun A character string specifying the covariance function.
//' Currently only works with "ard".
//' @param delta Scalar value that functions as a fixed nugget used to
//' stabilize matrix inverses and decompositions.
//' @param lnames Vector of the names of the length scale parameters, i.e.
//' "l1", ..., "ld"
//' @return Covariance matrix at the observed data locations if
//' x_pred = matrix(), or the covariances between x and x_pred if
//' x_pred is not an empty matrix.
//' @useDynLib sparseRGPs
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericMatrix make_cov_mat_ardC(NumericMatrix x,
                            NumericMatrix x_pred,
                            List cov_par,
                            String cov_fun,
                            double delta,
                            StringVector lnames){
  // select the squared exponential covariance function
  if(cov_fun == "ard"){

    // if the length of x_pred is zero, make variance covariance matrix for x
    if(NumericMatrix::is_na(x_pred(0,0)) == 1){
      int nrow = x.nrow();
      int ncol = x.nrow();

      NumericMatrix mat(nrow, ncol);

      for(int i = 0; i < nrow; i++){

        for(int j = 0; j < ncol; j++){
          if(i == j){
            mat(i , j) = cov_fun_sqrd_exp_ardC(x(i, _), x(j, _),
                cov_par,
                lnames = lnames) + pow(as<double>(cov_par["tau"]),2) + delta;
          }
          else{
            mat(i , j) = cov_fun_sqrd_exp_ardC(x(i, _), x(j, _),
                cov_par,
                lnames = lnames);
          }
        }

      }
      return mat;
    }
    // if the length of x_pred is nonzero, make the covariance matrix where x is rows, x_pred is columns
    if(NumericMatrix::is_na(x_pred(0,0)) != 1){
      int nrow = x.nrow();
      int ncol = x_pred.nrow();

      NumericMatrix mat(nrow, ncol);

      for(int i = 0; i < nrow; i++){

        for(int j = 0; j < ncol; j++){
          mat(i , j) = cov_fun_sqrd_exp_ardC(x(i, _),x_pred(j,_),
              cov_par,  lnames = lnames);
        }

      }
      return mat;
    }
  }

  else{
    NumericMatrix mat(0,0);
    Rcerr << "Error: invalid covariance function";
    return mat;
  }
  NumericMatrix mat(0,0);
  Rcerr << "Error: invalid covariance function";
  return mat;
}

/*** R
# x1 <- matrix(nrow = 2000, ncol = 2, rnorm(n = 12000))
# x2 <- matrix(nrow = 100, ncol = 2, rnorm(n = 200))
# print(x1)
# print(x2)
# cov_par <- list("sigma" = 1, "l" = 1, "tau" = 0)

# make_cov_matC(x = x1, x_pred = x2, cov_par = cov_par, cov_fun = "sqexp")

# system.time(test_cov_matC <- make_cov_matC(x = x1, x_pred = matrix(), cov_par = cov_par, cov_fun = "sqexp"))
#
# make_cov_matC(x = x1[1:10,], x_pred = matrix(), cov_par = cov_par, cov_fun = "exp")
*/
