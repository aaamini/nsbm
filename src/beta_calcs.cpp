// [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <random>
// #include <Rcpp/Benchmark/Timer.h>

// #include <boost/math/special_functions/factorials.hpp>
// #include <boost/format.hpp> 

// #include "sampling.h" 
// #include "utils.h" 

using namespace Rcpp;

// [[Rcpp::export]]
double comp_beta_ratio_v1(
    const double alpha, const double beta, 
    const int d, const int dbar) {

    return R::beta(alpha + d, beta + dbar) /
            // (R::beta(alpha, beta) + DBL_MIN);
             (R::beta(alpha, beta) + DBL_MIN);
}

// [[Rcpp::export]]
void print_smallest_double() {
    Rcout << "minimum double = " << DBL_MIN << "\n";
}

// // [[Rcpp::export]]
// double comp_beta_ratio_v2(
//     const double alpha, const double beta, 
//     const int d, const int dbar) {

//     // This still overflows    
//     return boost::math::rising_factorial(alpha, d) * 
//         boost::math::rising_factorial(beta, dbar) /
//         boost::math::rising_factorial(alpha + beta, d + dbar);
    
// }

// [[Rcpp::export]]
double comp_log_gamma_ratio(const double x, const int d) {
    // computes log [Gamma(x + d) / Gamma(x)]
    if (d == 0) return 0;
    if (d < 0) return -comp_log_gamma_ratio(x+d, -d);
    if (x <= 0) throw std::invalid_argument(
        "Invalid x value (<= 0) for Gamma ratio: " + std::to_string(x));
    // if (x <= 0) throw std::invalid_argument(str(boost::format("Invalid x value of %1% for Gamma ratio.") % x));

    // This is a more numerically stable way to compute the log gamma ratio
    double result = 0;
    for (double i = x; i < x+d; i++) {
        result += log(i);
    }
    return result;
}

// [[Rcpp::export]]
double comp_log_beta_ratio(
    const double alpha, const double beta, 
    const int d, const int dbar) {

    if (alpha <= 0) throw std::invalid_argument(
        "Invalid alpha value of " + std::to_string(alpha) + " for Beta ratio.");
    if (beta <= 0) throw std::invalid_argument(
        "Invalid beta value of " + std::to_string(beta) + " for Beta ratio.");
    if (alpha + d <= 0) throw std::invalid_argument(
        "Invalid alpha + d value of " + std::to_string(alpha + d) + " for Beta ratio.");
    if (beta + dbar <= 0) throw std::invalid_argument(
        "Invalid beta + dbar value of " + std::to_string(beta + dbar) + " for Beta ratio.");
        
    return comp_log_gamma_ratio(alpha, d) + 
        comp_log_gamma_ratio(beta, dbar) -
        comp_log_gamma_ratio(alpha + beta, d + dbar);
}