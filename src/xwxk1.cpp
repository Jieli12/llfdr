#include <Rcpp.h>
using namespace Rcpp;
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

//' @title Matrix Multiplication
//'
//' @param x matrix
//' @param w vector
//' @return matrix
//' @export
// [[Rcpp::export]]
NumericMatrix xwxk1(NumericMatrix x, NumericVector w) {
    int nrow = x.nrow(), ncol = x.ncol();
    NumericMatrix out(nrow, nrow);

    for (int i = 0; i < nrow; i++) {
        for (int k = 0; k < ncol; k++) {
            // double total = 0;
            double wx = 0;
            wx = w[k] * x(i, k);
            int j;
            for (j = 0; j < nrow; j++) {
                out(i, j) += wx * x(j, k);
            }
            // out(i, j) *= wx;
        }
    }
    return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
//
// /*** R
// xwxk1(Y, w)
// */
