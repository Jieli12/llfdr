#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//' @title LMEBeta2C
//' @description This function is used to compute the delete-i beta goal
//' function f = (sbuk / sbk - suk / sk)^2;
//'
//' @param v the beta when the i-th observation deleted
//' @param y2 the square of y
//' @param expu the exp(u_j - u_i)
//' @param k the kernel
//' @param uk the (u_j - u_i) * kernel
//' @param SK summation of kernel
//' @param SUK summation of kernel * u
//' @export
//' @return goal function value
// [[Rcpp::export]]
double LMEBeta2C(double v, NumericVector y2, NumericVector expu,
                NumericVector k, NumericVector uk, double SK, double SUK) {
    // int n = y2.length();
    // double temp1 = 0;
    // double ssbk = 0;
    // double ssbuk = 0;
    NumericVector temp1 = y2 / pow(expu, v);
    double sbk = sum(temp1 * k);
    double sbuk = sum(temp1 * uk);
    // for (int i = 0; i < n; i++) {
    //     ssbk += sbk[i];
    //     ssbuk += sbuk[i];
    // }
    return(pow(sbuk / sbk - SUK / SK, 2.0));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
//
// /*** R
// timesTwo(42)
// */
