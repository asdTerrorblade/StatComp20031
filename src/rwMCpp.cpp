#include <Rcpp.h>
using namespace Rcpp;

//' @title rwMCpp
//' @description A Metropolis sampler using Rcpp
//' @param N the number of samples
//' @param sigma variance
//' @param x0 the initial value
//' @return random walk samples
//' @examples
//' \dontrun{
//' rnC <- rwMCpp(1, 20, 1000)
//' plot(rnC)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMCpp (double sigma, double x0, int N) {
    NumericVector x(N);
    x[0] = x0;
    NumericVector u = runif(N);
    for (int i = 1; i < N;i++ ) {
        NumericVector y = rnorm(1, x[i-1], sigma);
        if (u[i] <= (exp(-abs(y[0])) / exp(-abs(x[i-1])))){
            x[i] = y[0];
        }
                else {
            x[i] = x[i-1];
        }
    }
    return(x);
}