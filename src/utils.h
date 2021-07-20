
#ifndef __UTILS__
#define __UTILS__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

int find_tunc(arma::vec beta, double threshold);

arma::uvec get_freq(arma::uvec z, int K);
arma::uvec get_up_freq(arma::uvec freq);

arma::umat get_freq_minus_self(arma::uvec z, int K);

arma::vec fast_agg(arma::vec x, arma::uvec z, int K);

arma::uvec fast_agg_u(arma::uvec x, arma::uvec z, int K) ;

arma::mat comp_blk_sums(arma::sp_mat At, arma::uvec z, int Kcap);

arma::mat sp_compress_col(arma::sp_mat At, arma::uvec z, int Kcap);

List comp_blk_sums_and_sizes(arma::sp_mat At, arma::uvec z, int Kcap, bool div_diag = true);

arma::vec sp_single_col_compress(arma::sp_mat A, int col_idx, arma::uvec z, int Kcap);

arma::mat comp_blk_sums_diff(arma::sp_mat& A, int s, int zs_new, arma::uvec& z, int Kcap);

arma::mat beta_fun_symmat(arma::mat a, arma::mat b);

arma::mat comp_beta_matrix(const arma::sp_mat& A, arma::uvec& z, const int K, double alpha, double beta);

arma::mat comp_blk_sums_diff_v1(const arma::vec& U, const int zs_new, const int zs_old);
arma::mat comp_blk_sums_diff_v2(const arma::vec& U, const int zs_new, const int zs_old);

arma::vec comp_beta_ratio_prods_v1(
    const arma::mat& m, 
    const arma::mat& mbar, 
    const arma::vec& U,
    const arma::uvec& V, 
    const int zs_old,
    const int alpha, const int beta);

arma::vec comp_log_beta_ratio_sums(
    const arma::mat& m, 
    const arma::mat& mbar, 
    const arma::vec& U,
    const arma::uvec& V, 
    const int zs_old,
    const int alpha, const int beta);

    
void print_progress(int itr, int itr_max);                    
#endif /* __UTILS__ */
