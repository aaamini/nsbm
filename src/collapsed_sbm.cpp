// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

#include "sampling.h" 
#include "utils.h" 

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec f_fun(arma::vec x, double alpha, double beta) {
    int n =  x.n_elem;
    arma::vec out = arma::vec(n, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
        out(i) = ::tgamma(x(i) + alpha) / exp(log(beta+0.5)*(x(i) + alpha));
    }
    return out;
    // return Rcpp::as<arma::vec>( 
    //     Rcpp::gamma(Rcpp::wrap(x + alpha)) / (exp(log(beta+0.5)*(Rcpp::wrap(x + alpha))
    // );
}

// [[Rcpp::export]]
double ratio_fun(arma::vec N, arma::vec Np, double alpha, double beta) { // , int r, int rp) {
    // double out = 0;
    // int K = N.n_rows;
    return arma::prod(f_fun(Np/2, alpha, beta) /
                        f_fun(N/2, alpha, beta));
}


//' @export
// [[Rcpp::export]]
List fit_dcsbm(arma::sp_mat& A, const int K, 
                    const double alpha = 1, const double beta = 1,
                    const int niter = 2) {

    int n = A.n_rows;
    // arma::umat z_hist(n, niter*n, arma::fill::zeros);
    arma::umat z_hist(n, niter, arma::fill::zeros);
    arma::umat nn_hist(K, niter, arma::fill::zeros);
    arma::mat  N_hist(K*K, niter, arma::fill::zeros);
     

    // Randomly initialize labels
    arma::uvec z = sample_int_vec(K, n);
    arma::uvec nn = get_freq(z, K);
    arma::mat  N = comp_blk_sums(A, z, K);
    
    for (int iter = 0; iter < niter; iter++) {
        z_hist.col(iter) = z;
        nn_hist.col(iter) = nn;
        N_hist.col(iter) = N.as_col();

        for (int s = 0; s < n; s++) {
            // z_hist.col(iter*n+s) = z;
            arma::vec prob(K, arma::fill::zeros);
            arma::mat Np(K, K, arma::fill::zeros);

            for (int rp = 0; rp < K; rp++) { // rp is the potential new value of z(s)
                Np = N + comp_blk_sums_diff(A, s, rp, z, K);
                // Rcpp::Rcout <<  "s = " << s << "\n" << comp_blk_sums_diff(A, s, rp, z, K) << "\n";
                prob(rp) = ratio_fun(N.as_col(), Np.as_col(), alpha, beta);
                prob(rp) *= static_cast<double>((nn(rp) + 1)) / nn(z(s)); 
            } // rp
          

            int old_zs = z(s); // previous value of z(s)
            int new_zs = sample_index(prob);
            //Rcpp::Rcout <<  "s = " << s << "\n" << z.t();

            // Rcpp::Rcout << prob /sum(prob) << "\n" << "old = " << old_zs << "new = " << z(s) << "\n";

            // z is still old z    
            N = N + comp_blk_sums_diff(A, s, new_zs, z, K); // update N
            z(s) = new_zs; // update z
            //Rcpp::Rcout <<  "s = " << s << "\n" << comp_blk_sums_diff(A, s, z(s), z, K) << "\n";
            // update nn
            nn(old_zs)--;
            nn(new_zs)++;
        } // s
    }
    return Rcpp::List::create( 
        Rcpp::Named("z_hist") = z_hist,
        Rcpp::Named("nn_hist") = nn_hist,
        Rcpp::Named("N_hist") = N_hist
    ); 
}
