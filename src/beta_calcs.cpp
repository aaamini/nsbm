// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
// #include <Rcpp/Benchmark/Timer.h>

// #include "sampling.h" 
// #include "utils.h" 

using namespace Rcpp;

// [[Rcpp::export]]
double comp_beta_ratio_v1(
    const double alpha, const double beta, 
    const int d, const int dbar) {

    return R::beta(alpha + d, beta + dbar) /
            (R::beta(alpha, beta) + DBL_MIN);
}

// [[Rcpp::export]]
void print_smallest_double() {
    Rcout << "minimum double = " << DBL_MIN << "\n";
}
