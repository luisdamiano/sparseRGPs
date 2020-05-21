#include <Rcpp.h>
//#include <boost/lexical_cast.hpp>
using namespace Rcpp;


// inverse transformation functions
//' @export
// [[Rcpp::export]]
NumericVector real_to_pos(NumericVector x){
  return exp(x);
}

//' @export
// [[Rcpp::export]]
NumericVector pos_to_real(NumericVector x){
  return log(x);
}

//' @export
// [[Rcpp::export]]
NumericVector real_to_bounded(NumericVector x, NumericVector ub, NumericVector lb){
  return (ub * exp(x) + lb) / (exp(x) + 1);
}


// Squared exponential derivatives

// [[Rcpp::export]]
Rcpp::List dsqexp_dsigmaC(NumericVector x1, NumericVector x2, List cov_par){

  double sigma = as<double>(cov_par["sigma"]);
  double l = as<double>(cov_par["l"]);
  double dsigma_dsigmat = sigma;

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");
  NumericVector sigma_vec(1);
  sigma_vec[0] = sigma;

  return Rcpp::List::create(
    Rcpp::Named("derivative") = 2 * sigma * exp( -(1/(2*pow(l,2))) * sum(pow((x1 - x2), 2) )) * dsigma_dsigmat,
    Rcpp::Named("trans_par") = log(sigma),
    Rcpp::Named("inv_trans_par") = real_to_pos(sigma)
  );

}

// [[Rcpp::export]]
Rcpp::List dsqexp_dsigma_ardC(NumericVector x1,
                              NumericVector x2,
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

  double dsigma_dsigmat = sigma;

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");
  NumericVector sigma_vec(1);
  sigma_vec[0] = sigma;

  return Rcpp::List::create(
    Rcpp::Named("derivative") = 2 * sigma * exp( -(sum(pow((x1 - x2) / l, 2) ) / 2) ) * dsigma_dsigmat,
    Rcpp::Named("trans_par") = log(sigma),
    Rcpp::Named("inv_trans_par") = real_to_pos(sigma)
  );

}

// [[Rcpp::export]]
Rcpp::List dsqexp_dlC(NumericVector x1, NumericVector x2, List cov_par){

  double sigma = as<double>(cov_par["sigma"]);
  double l = as<double>(cov_par["l"]);
  double dl_dlt = l;

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");
  NumericVector l_vec(1);
  l_vec[0] = l;

  return Rcpp::List::create(
    Rcpp::Named("derivative") = (pow(sigma,2) * exp( (-1/(2*pow(l,2))) *
      sum( pow((x1 - x2),2)) )) * ( (1/(pow(l,3))) * sum(pow((x1 - x2),2) )) * dl_dlt,
    Rcpp::Named("trans_par") = log(l),
    Rcpp::Named("inv_trans_par") = real_to_pos(l)
  );

}

// [[Rcpp::export]]
Rcpp::List dsqexp_dl_ardC(NumericVector x1, NumericVector x2,
                          List cov_par,
                          StringVector lnames,
                          double comp){

  double sigma = as<double>(cov_par["sigma"]);
  int d = x1.size();
  NumericVector l(d);

  for(int i = 0; i < d; i++)
  {
    String temp = lnames[i];
    l[i] = as<double>(cov_par[temp]);
  }
  comp = comp - 1;
  double dl_dlt = l[comp];
  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");

  //Rcout << l << "\n";
  //Rcout << real_to_pos(l) << "\n";

  //NumericVector l_vec(1);
  //l_vec[0] = l[comp];

  return Rcpp::List::create(
    Rcpp::Named("derivative") = (pow(sigma,2) * exp( -(sum(pow((x1 - x2) / l, 2) ) / 2) ) ) *
      ( (1/(pow(l[comp],3))) * (pow((x1[comp] - x2[comp]),2) )) * dl_dlt,
      Rcpp::Named("trans_par") = log(l[comp]),
      Rcpp::Named("inv_trans_par") = real_to_pos(l[comp])
  );

}

// [[Rcpp::export]]
Rcpp::List dsqexp_dtauC(NumericVector x1, NumericVector x2, List cov_par){

  //double sigma = as<double>(cov_par["sigma"]);
  //double l = as<double>(cov_par["l"]);
  double tau = as<double>(cov_par["tau"]);
  double dtau_dtaut = tau;

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");
  NumericVector tau_vec(1);
  tau_vec[0] = tau;

  //LogicalVector check = x1 == x2;

  double deriv;
  if(is_true(all(x1 == x2)))
  {
    deriv = 2 * tau * dtau_dtaut;
  }
  else{
    deriv = 0;
  }

  return Rcpp::List::create(
    Rcpp::Named("derivative") = deriv,
    Rcpp::Named("trans_par") = log(tau),
    Rcpp::Named("inv_trans_par") = real_to_pos(tau_vec)

  );
}



// [[Rcpp::export]]
Rcpp::List dsqexp_dx2C(NumericVector x1, NumericVector x2, List cov_par, NumericVector lb, NumericVector ub){

  double sigma = as<double>(cov_par["sigma"]);
  double l = as<double>(cov_par["l"]);
  NumericVector tx2 = log((x2 - lb) / (ub - x2));
  NumericVector dx2_dtx2 = (exp(tx2) * (ub - lb)) / pow((exp(tx2) + 1), 2);

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");


  return Rcpp::List::create(
    Rcpp::Named("derivative") = (1/(pow(l,2))) * (x1 - x2) * pow(sigma,2) *
      (exp( -sum( pow((x1 - x2),2)) / (2*pow(l,2)))) * dx2_dtx2,
    Rcpp::Named("trans_par") = tx2,
    Rcpp::Named("inv_trans_par") = real_to_bounded(x2, ub = ub, lb = lb)
  );

}

// [[Rcpp::export]]
Rcpp::List dsqexp_dx2_ardC(NumericVector x1, NumericVector x2,
                           List cov_par,
                           NumericVector lb,
                           NumericVector ub,
                           StringVector lnames){

  double sigma = as<double>(cov_par["sigma"]);
  int d = x1.size();
  NumericVector l(d);

  for(int i = 0; i < d; i++)
  {
    String temp = lnames[i];
    l[i] = as<double>(cov_par[temp]);
  }

  NumericVector tx2 = log((x2 - lb) / (ub - x2));
  NumericVector dx2_dtx2 = (exp(tx2) * (ub - lb)) / pow((exp(tx2) + 1), 2);

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");


  return Rcpp::List::create(
    Rcpp::Named("derivative") = (1/(pow(l,2))) * (x1 - x2) * pow(sigma,2) *
      (exp( -(sum(pow((x1 - x2) / l, 2) ) / 2) )) * dx2_dtx2,
      Rcpp::Named("trans_par") = tx2,
      Rcpp::Named("inv_trans_par") = real_to_bounded(x2, ub = ub, lb = lb)
  );

}

// exponential covariance function derivatives

// [[Rcpp::export]]
Rcpp::List dexp_dsigmaC(NumericVector x1, NumericVector x2, List cov_par){

  double sigma = as<double>(cov_par["sigma"]);
  double l = as<double>(cov_par["l"]);
  double dsigma_dsigmat = sigma;

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");
  NumericVector sigma_vec(1);
  sigma_vec[0] = sigma;

  return Rcpp::List::create(
    Rcpp::Named("derivative") = 2 * sigma * exp( -(1/(l)) * sqrt(sum( pow((x1 - x2),2) ))) * dsigma_dsigmat,
    Rcpp::Named("trans_par") = log(sigma),
    Rcpp::Named("inv_trans_par") = real_to_pos(sigma)
  );

}

// [[Rcpp::export]]
Rcpp::List dexp_dlC(NumericVector x1, NumericVector x2, List cov_par){

  double sigma = as<double>(cov_par["sigma"]);
  double l = as<double>(cov_par["l"]);
  double dl_dlt = l;

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");
  NumericVector l_vec(1);
  l_vec[0] = l;

  return Rcpp::List::create(
    Rcpp::Named("derivative") = (pow(sigma,2) * exp( (-1/(l)) * sqrt(sum( pow((x1 - x2),2) )) )) * ( (1/(pow(l,2))) * sqrt(sum( pow((x1 - x2),2) )) ) * dl_dlt,
    Rcpp::Named("trans_par") = log(l),
    Rcpp::Named("inv_trans_par") = real_to_pos(l)
  );

}

// [[Rcpp::export]]
Rcpp::List dexp_dtauC(NumericVector x1, NumericVector x2, List cov_par){

  //double sigma = as<double>(cov_par["sigma"]);
  //double l = as<double>(cov_par["l"]);
  double tau = as<double>(cov_par["tau"]);
  double dtau_dtaut = tau;

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");
  NumericVector tau_vec(1);
  tau_vec[0] = tau;

  double deriv;
  if(is_true(all(x1 == x2)))
  {
    deriv = 2 * tau * dtau_dtaut;
  }
  else{
    deriv = 0;
  }

  return Rcpp::List::create(
    Rcpp::Named("derivative") = deriv,
    Rcpp::Named("trans_par") = log(tau),
    Rcpp::Named("inv_trans_par") = real_to_pos(tau)

  );


}


//matrix population function

// [[Rcpp::export]]
NumericMatrix dsig_dthetaC(NumericMatrix x, NumericMatrix x_pred,
                           List cov_par, String cov_fun, String par_name){
  // all positive only variables are transformed to the whole real line
  // and the gradient of the transformed variable + the untransformed current value
  // are returned

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");

  // select the squared exponential covariance function
  if(cov_fun == "sqexp"){
    int nrow = x.nrow();
    int ncol = x.nrow();

    // if the length of x_pred is zero, take derivatives of the covariance matrix for x
    if(NumericMatrix::is_na(x_pred(0,0)) == 1){

      NumericMatrix mat(nrow, ncol);

      // choose the parameter - sigma
      if(par_name == "sigma"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dsqexp_dsigmaC(x(i, _), x(j, _), cov_par = cov_par)["derivative"];

          }

        }
        return mat;
      }

      // choose the parameter - l
      if(par_name == "l"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dsqexp_dlC(x(i, _), x(j, _), cov_par = cov_par)["derivative"];

          }

        }
        return mat;
      }

      // choose the parameter - tau
      if(par_name == "tau"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dsqexp_dtauC(x(i, _), x(j, _), cov_par = cov_par)["derivative"];

          }

        }
        return mat;
      }

      //return mat;
    }
    // if the length of x_pred is nonzero, make the covariance matrix where x is rows, x_pred is columns
    if(NumericMatrix::is_na(x_pred(0,0)) != 1){
      int nrow = x.nrow();
      int ncol = x_pred.nrow();

      NumericMatrix mat(nrow, ncol);

      //choose covariance parameter - sigma
      if(par_name == "sigma")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dsqexp_dsigmaC(x(i, _), x_pred(j, _), cov_par = cov_par)["derivative"];
          }

        }
        return mat;
      }

      //choose covariance parameter - l
      if(par_name == "l")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dsqexp_dlC(x(i, _), x_pred(j, _), cov_par = cov_par)["derivative"];
          }

        }
        return mat;
      }

      //choose covariance parameter - tau
      if(par_name == "tau")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dsqexp_dtauC(x(i, _), x_pred(j, _), cov_par = cov_par)["derivative"];
          }

        }
        return mat;
      }
      else{
        NumericMatrix mat(0,0);
        Rcerr << "Error: invalid parameter name for chosen covariance function";
        return mat;
      }

      return mat;
    }

  }



  // select the exponential covariance function
  if(cov_fun == "exp"){

    // if the length of x_pred is zero, take derivatives of the covariance matrix for x
    if(NumericMatrix::is_na(x_pred(0,0)) == 1){
      int nrow = x.nrow();
      int ncol = x.nrow();

      NumericMatrix mat(nrow, ncol);

      // choose the parameter - sigma
      if(par_name == "sigma"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dexp_dsigmaC(x(i, _), x(j, _), cov_par = cov_par)["derivative"];

          }

        }
        return mat;
      }

      // choose the parameter - l
      if(par_name == "l"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dexp_dlC(x(i, _), x(j, _), cov_par = cov_par)["derivative"];

          }

        }
        return mat;
      }

      // choose the parameter - tau
      if(par_name == "tau"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dexp_dtauC(x(i, _), x(j, _), cov_par = cov_par)["derivative"];

          }

        }
        return mat;
      }
    }
    // if the length of x_pred is nonzero, make the covariance matrix where x is rows, x_pred is columns
    if(NumericMatrix::is_na(x_pred(0,0)) != 1){
      int nrow = x.nrow();
      int ncol = x_pred.nrow();

      NumericMatrix mat(nrow, ncol);

      //choose covariance parameter - sigma
      if(par_name == "sigma")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dexp_dsigmaC(x(i, _), x_pred(j, _), cov_par = cov_par)["derivative"];
          }

        }
        return mat;

      }

      //choose covariance parameter - l
      if(par_name == "l")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dexp_dlC(x(i, _), x_pred(j, _), cov_par = cov_par)["derivative"];
          }

        }
        return mat;
      }

      return mat;

      //choose covariance parameter - tau
      if(par_name == "tau")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dexp_dtauC(x(i, _), x_pred(j, _), cov_par = cov_par)["derivative"];
          }

        }
        return mat;
      }

      else{
        NumericMatrix mat(0,0);
        Rcerr << "Error: invalid parameter name for chosen covariance function";
        return mat;
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
  Rcerr << "Error";
  return mat;
}

// [[Rcpp::export]]
NumericMatrix dsig_dtheta_ardC(NumericMatrix x,
                               NumericMatrix x_pred,
                               List cov_par,
                               String cov_fun,
                               String par_name,
                               StringVector lnames){
  // all positive only variables are transformed to the whole real line
  // and the gradient of the transformed variable + the untransformed current value
  // are returned

  Function real_to_pos("real_to_pos");
  Function real_to_bounded("real_to_bounded");

  // select the squared exponential covariance function
  if(cov_fun == "ard"){
    int nrow = x.nrow();
    int ncol = x.nrow();

    // if the length of x_pred is zero, take derivatives of the covariance matrix for x
    if(NumericMatrix::is_na(x_pred(0,0)) == 1){

      NumericMatrix mat(nrow, ncol);

      // choose the parameter - sigma
      if(par_name == "sigma"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dsqexp_dsigma_ardC(x(i, _), x(j, _),
                cov_par = cov_par,
                lnames = lnames)["derivative"];

          }

        }
        return mat;
      }

      // choose the parameter - l
      LogicalVector lname_test(lnames.size());
      int comp;
      for(int i = 0; i < lnames.size(); i++)
      {
        lname_test[i] = (lnames[i] == par_name);
        if(lname_test[i] == TRUE)
        {
          comp = i + 1;
        }
      }

      //bool jim = any(lname_test == TRUE);
      //Rcout << "comp is " << comp << "\n";
      //Rcout << jim << "\n";
      if(is_true(any(lname_test == TRUE))){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dsqexp_dl_ardC(x(i, _), x(j, _), cov_par = cov_par,
                lnames = lnames, comp = comp)["derivative"];

          }

        }


        return mat;
      }

      // choose the parameter - tau
      if(par_name == "tau"){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dsqexp_dtauC(x(i, _), x(j, _), cov_par = cov_par)["derivative"];

          }

        }
        return mat;
      }

      //return mat;
    }
    // if the length of x_pred is nonzero, make the covariance matrix where x is rows, x_pred is columns
    if(NumericMatrix::is_na(x_pred(0,0)) != 1){
      int nrow = x.nrow();
      int ncol = x_pred.nrow();

      NumericMatrix mat(nrow, ncol);

      //choose covariance parameter - sigma
      if(par_name == "sigma")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dsqexp_dsigma_ardC(x(i, _), x_pred(j, _),
                cov_par = cov_par, lnames = lnames)["derivative"];
          }

        }
        return mat;
      }

      // choose the parameter - l
      LogicalVector lname_test(lnames.size());
      int comp;
      for(int i = 0; i < lnames.size(); i++)
      {
        lname_test[i] = (lnames[i] == par_name);
        if(lname_test[i] == TRUE)
        {
          comp = i + 1;
        }
      }

      if(is_true(any(lname_test == TRUE))){

        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){

            mat(i , j) = dsqexp_dl_ardC(x(i, _), x_pred(j, _), cov_par = cov_par,
                lnames = lnames, comp = comp)["derivative"];

          }

        }
        return mat;
      }

      //choose covariance parameter - tau
      if(par_name == "tau")
      {
        for(int i = 0; i < nrow; i++){

          for(int j = 0; j < ncol; j++){
            mat(i , j) = dsqexp_dtauC(x(i, _), x_pred(j, _), cov_par = cov_par)["derivative"];
          }

        }
        return mat;
      }
      else{
        NumericMatrix mat(0,0);
        Rcerr << "Error: invalid parameter name for chosen covariance function";
        return mat;
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
  Rcerr << "Error";
  return mat;
}

/*** R
# x1 <- matrix(nrow = 5, ncol = 2, rnorm(n = 10))
# x2 <- matrix(nrow = 1, ncol = 2, rnorm(n = 2))
# # print(x1)
# # print(x2)
# cov_par <- list("sigma" = 1, "l" = 1, "tau" = 0)
#
# dsig_dthetaC(x = x1, x_pred = matrix(), cov_par = cov_par, cov_fun = "sqexp", par_name = "l")
# dsig_dthetaC(x = x1, x_pred = x2, cov_par = cov_par, cov_fun = "sqexp", par_name = "sigma")
# dsig_dthetaC(x = x1, x_pred = x2, cov_par = cov_par, cov_fun = "sqexp", par_name = "l")
# dsig_dthetaC(x = x1, x_pred = x2, cov_par = cov_par, cov_fun = "exp", par_name = "l")

*/
