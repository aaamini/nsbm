// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "utils.h"
#include "sampling.h"
#include "beta_calcs.h"
#include "dpsbm.h"

using namespace Rcpp;

struct BetaParameters {
  double alpha;
  double beta;
  
  BetaParameters() :alpha{1}, beta{1} {};
  BetaParameters(const double a, const double b) :alpha{a}, beta{b} {};
};


// template<typename T>
// T vector_lbeta(T A, T B) {
//     for(auto& x : vec_obj) {
//         x = R::lbeta(x);
//     }
//     return vec_obj;
// }

// [[Rcpp::export]]
arma::cube cube_lbeta(const arma::cube& A, const arma::cube& B) {
  arma::cube result(arma::size(A));
  for (int i = 0; i < A.n_rows; i++) {
    for (int j = 0; j < A.n_cols; j++) {
      for (int k = 0; k < A.n_cols; k++) {
        result(i,j,k) = R::lbeta(A(i,j,k), B(i,j,k));
      }
    }
  }
  return result;
}

// 
// arma::cube test_vector_lbeta(arma::cube X) {
//     return vector_lbeta(X);
// }

class NestedSBM {
public:
  List A;
  int J;  // number of networks
  int K;  // truncation level for the prior on "z"
  int L;  // truncation level for the prior on "w"
  arma::uvec n; // n(j) = number nodes in network j
  arma::uvec z; // z(j) = label of network j
  std::vector<arma::uvec> xi; // xi[j](s) = label of nodes s in network j
  arma::vec pi;   // the prior on "z"
  arma::mat w;    // the prior on "xi"
  arma::vec u;
  arma::mat v;
  // std::vector<arma::mat> blk_compressions;
  
  arma::cube a;   // log[eta/(1-eta)]
  arma::cube b;   // log(1-eta)
  arma::cube eta;
  
  arma::cube m;   // m(x,y,k)
  arma::cube mbar; // mbar(x,y,k)
  
  arma::uvec z_freq;  // In the current implementation, these two parameters only hold correct value after update_w and upate_pi
  arma::umat xi_freq_over_z; // TODO: make them hold correct valu after every label update
  
  BetaParameters beta_params;
  double w0;
  double pi0;
  
  // arma::vec z_log_prob_record; // for diagnostics
  
  NestedSBM(List A_, 
            // const arma::uvec z_init, 
            const int K, const int L, 
            const double alpha_eta = 1, 
            const double beta_eta = 1,
            const double w0_init = 0, 
            const double pi0_init = 0) : K{K}, L{L}, w0{w0_init}, pi0{pi0_init} {
              
              // initialize and allocate variables
              A = A_;
              J = A.length();
              n = arma::uvec(J);
              xi = std::vector<arma::uvec>(J);
              // blk_compressions = std::vector<arma::mat>(J);
              // z = arma::uvec(J);
              m = arma::cube(L, L, K, arma::fill::zeros);     // m tensor
              mbar = arma::cube(L, L, K, arma::fill::zeros);  // mbar tensor
              a = arma::cube(L, L, K, arma::fill::zeros);     // a tensor
              b = arma::cube(L, L, K, arma::fill::zeros);     // b tensor
              eta = arma::cube(L, L, K, arma::fill::zeros); 
              // pi = arma::vec(K, arma::fill::ones);
              w = arma::mat(L, K, arma::fill::ones);
              
              // z_log_prob_record = arma::vec(K);
              
              for (int j=0; j < J; j++) {
                n(j) = Rcpp::as<arma::sp_mat>(A[j]).n_rows;
                // xi[j] = sample_int_vec(L, n[j]);
              }
              // z = sample_int_vec(K, J);
              beta_params.alpha = alpha_eta;
              beta_params.beta = beta_eta;
              set_xi_to_random_labels();
              set_z_to_random_labels();
              
              z_freq = get_freq(z, K);
              xi_freq_over_z = arma::umat(L, K, arma::fill::ones);
              
              // Rcout << "1";
              for (int k = 0; k < K; k++) {
                xi_freq_over_z.col(k) = get_xi_freq_over_z(k);
              }
              
              // Rcout << "2";
              if (pi0 == 0) pi0 = 1. / (J * log(J));
              if (w0 == 0)  {
                int n_min = arma::min(n);
                w0 = 1. / (n_min * log(n_min));
              }
              
              // Rcout << "3";
              arma::vec u = arma::vec(K, arma::fill::ones);
              //u(arma::span(0,K-2)) *= 1. / (1+pi0);            
              u(arma::span(0,K-2)) *= 0.5;            
              pi = stick_break(u);
              // pi = arma::vec(K, arma::fill::ones);
              
              // Rcout << "4";
              v = arma::mat(L, K, arma::fill::ones);
              w = arma::mat(L, K, arma::fill::ones);
              v.rows(0,L-2) *= 1. /(1+w0);
              // v.rows(0,K-2) *= 0.5;
              // Rcpp::print(wrap(v));
              
              // Rcout << "5";
              for (int k = 0; k < K; k++) {
                w.col(k) = stick_break(v.col(k));
              }
              // w = arma::mat(L, K, arma::fill::ones);
              
              // for (int j = 0; j < J; j++) {
              //     blk_compressions[j] = sp_compress_col(A[j], xi[j], L);
              // }
              // Rcpp::print(wrap( blk_compressions[1]));
            }
  
  // --- updates for non-collapsed sampler --->
  void update_xi_element_via_eta(const int j, const int s) {
    arma::vec tau = sp_single_col_compress(A[j], s, xi[j], L);
    arma::uvec rho = get_freq(xi[j], L);
    rho(xi[j](s))--;         
    
    int xi_j_s_old = xi[j](s);
    arma::vec log_prob = log(w.col(z(j))) + a.slice(z(j)) * tau + b.slice(z(j)) * rho; 
    xi[j](s) = sample_index(safe_exp(log_prob));
    
    // update m and mbar
    //if ( xi_j_s_old != xi[j](s) ) {
    //  List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
    //  arma::mat lambda = out["lambda"];
    //  arma::umat NN = out["NN"]; 
    //  m.slice(z(j)) += lambda;
    //  mbar.slice(z(j)) += NN - lambda;
    //}
    
  }
  
  void update_xi_element_via_marginal(const int j, const int s) {
    arma::vec tau = sp_single_col_compress(A[j], s, xi[j], L);
    arma::uvec rho = get_freq(xi[j], L);
    rho(xi[j](s))--;
    
    int xi_j_s_old = xi[j](s);
    arma::vec prob(L, arma::fill::zeros);
    for (int l = 0; l < L; l++) {
      for (int k = 0; k < K; k++) {
        prob(l) += w.col(k)(l) * safe_exp(a.slice(k) * tau + b.slice(k) * rho)(l) * pi(k);
      }
    }
    
    xi[j](s) = sample_index(prob);
    
    // update m and mbar
    //if ( xi_j_s_old != xi[j](s) ) {
    //  List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
    //  arma::mat lambda = out["lambda"];
    //  arma::umat NN = out["NN"]; 
    //  m.slice(z(j)) += lambda;
    //  mbar.slice(z(j)) += NN - lambda;
    //}
    
  }
  
  void comp_count_tensors() {
    //List out = comp_blk_sums_and_sizes(Rcpp::as<arma::sp_mat>(A[0]), xi[0], L);
    m = arma::cube(L, L, K, arma::fill::zeros); // m tensor
    mbar = arma::cube(L, L, K, arma::fill::zeros); // mbar tensor
    
    for (int j=0; j < J; j++) {
      List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
      arma::mat lambda = out["lambda"];
      arma::umat NN = out["NN"]; 
      m.slice(z(j)) += lambda;
      mbar.slice(z(j)) += NN - lambda;
    }
  }
  
  void update_eta() {    
    // Update the eta-related tensors
    comp_count_tensors();
    for (int k = 0; k < K; k++) {
      eta.slice(k) = symmat_rbeta(m.slice(k) + beta_params.alpha, mbar.slice(k) + beta_params.beta);
      a.slice(k) = log(eta.slice(k) / (1-eta.slice(k)) + perturb);
      b.slice(k) = log(1-eta.slice(k) + perturb);
    }
  }
  
  void update_eta_with_constraint() {    
    // Update the eta-related tensors
    comp_count_tensors();
    for (int k = 0; k < K; k++) {
      for (int x = 0; x < L; x++) {
        for (int y = x; y < L; y++) {
          
          double max_eta = 0;
          for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
              if (i == x & j == y) continue; 
              max_eta = std::max(eta(i,j,k), max_eta);
            }
          }
          
          //eta(x,y,k) = rTruncBeta_tail(m(x,y,k) + beta_params.alpha, mbar(x,y,k) + beta_params.beta, .5 * max_eta, 1);
          do {
            eta(x,y,k) = R::rbeta(m(x,y,k) + beta_params.alpha, mbar(x,y,k) + beta_params.beta);
          }
          while (eta(x,y,k) < .05 * max_eta);
          a(x,y,k) = log(eta(x,y,k) / (1-eta(x,y,k)) + perturb);
          b(x,y,k) = log(1-eta(x,y,k) + perturb);
          
          eta(y,x,k) = eta(x,y,k);
          a(y,x,k) = a(x,y,k);
          b(y,x,k) = b(x,y,k);
        }
      }
    }
  }
  
  void update_z_element_via_eta(const int j) {
    int r0 = z(j);
    List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
    arma::mat D = out["lambda"]; 
    arma::umat M = out["NN"];
    
    arma::uvec xi_j_freq = get_freq(xi[j], L);
    arma::vec log_prob = log(w.t() + perturb) * xi_j_freq;
    log_prob += log(pi + perturb);
    
    for (int r = 0; r < K; r++) {
      for (int x = 0; x < L; x++) {
        for (int y = x; y < L; y++) {
          log_prob(r) += 
            (a(x,y,r) - a(x,y,r0))*D(x,y) + 
            (b(x,y,r) - b(x,y,r0))*M(x,y);
        }
      }
    }
    
    // update z(j)
    z(j) = sample_index(safe_exp(log_prob));
    
    // update m and mbar tensors
    //if (z(j) != r0) {
      
      //comp_count_tensors();
      
      //m.slice(z(j)) += D;
      //mbar.slice(z(j)) += M - D;

      //m.slice(r0) -= D;
      //mbar.slice(r0) -= M - D;
    //}
    
  }
  
  void set_xi_to_dpsbm_labels(const int niter) {
    for (int j=0; j < J; j++) {
      xi[j] = fit_dpsbm(A[j], w0, beta_params.alpha, beta_params.beta, niter, L).col(niter-1)-1;
    }
  }
  
  // List run_gibbs(const int niter = 100, const bool init_count_tensors = true) {
  List run_gibbs_via_eta(const int niter, const int version) {
    // Run full Gibbs updates for "niter" iterations and record label history
    std::vector<std::vector<arma::uvec>> xi_hist(niter+1);
    arma::umat z_hist(J, niter+1);
    
    xi_hist[0] = xi;
    z_hist.col(0) = z + 1;
    
    
    
    // comp_count_tensors();
    for (int iter = 0; iter < niter; iter++) {
      
      if (iter == 0) comp_count_tensors();
      
      switch (version) {
      case 1:
        update_eta(); // also updates count tensors m and mbar 
        for (int j = 0; j < J; j++) {
          update_z_element_via_eta(j);
          for (int s = 0; s < n(j); s++) {
            update_xi_element_via_eta(j, s);
          } // s
        } // j
        
        update_pi();
        update_w();
        break;
        
      case 2:
        update_eta(); // also updates count tensors m and mbar 
        update_pi();
        update_w();
        
        for (int j = 0; j < J; j++) {
          update_z_element_via_eta(j);
          for (int s = 0; s < n(j); s++) {
            update_xi_element_via_eta(j, s);
          } // s
          
        } // j
        break;
        
      case 3:
        if (iter == 0) {
          set_xi_to_dpsbm_labels(50);
          xi_hist[0] = xi;
        }
        update_eta(); // also updates count tensors m and mbar 
        for (int j = 0; j < J; j++) {
          update_z_element_via_eta(j);
          for (int s = 0; s < n(j); s++) {
            update_xi_element_via_eta(j, s);
          } // s
        } // j
        
        update_pi();
        update_w();
        break;
        
      case 4:
        update_eta_with_constraint(); // also updates count tensors m and mbar 
        update_pi();
        update_w();
        
        for (int j = 0; j < J; j++) {
          update_z_element_via_eta(j);
          for (int s = 0; s < n(j); s++) {
            update_xi_element_via_eta(j, s);
          } // s
          
        } // j
        break;
        
      case 5:
        if (iter == 0) {
          set_xi_to_dpsbm_labels(50);
          xi_hist[0] = xi;
        }
        update_eta(); // also updates count tensors m and mbar
        for (int j = 0; j < J; j++) {
          for (int s = 0; s < n(j); s++) {
            update_xi_element_via_marginal(j, s);
          } // s
          update_z_element_via_eta(j);
        } // j
        
        update_pi();
        update_w();
        break;
        
      case 6:
        if (iter == 0) {
          set_xi_to_dpsbm_labels(50); 
          xi_hist[0] = xi;
        }
        update_eta(); // also updates count tensors m and mbar 
        for (int j = 0; j < J; j++) {
          update_z_element_via_eta(j);
          for (int s = 0; s < n(j); s++) {
            update_xi_element_via_marginal(j, s);
          } // s
        } // j
        
        update_pi();
        update_w();
        break;
        
      }        
      
      
      // Rcpp::print(wrap(z.t()));
      xi_hist[iter+1] = xi;
      z_hist.col(iter+1) = z + 1;
    } // iter
    
    return Rcpp::List::create( 
      Rcpp::Named("z") = z_hist,
      Rcpp::Named("xi") = xi_hist
    );
  }
  
  // <--- end of updates for non-collapsed sampler ---
  
  void print() {
    Rcout << "- Nested SMB Model - \n"
          << "  J = " << J 
          <<"\n  K = " << K 
          << "\n  L = " << L << "\n"
          << std::left 
          << std::setw(15) << "  (alpha, beta)" 
          <<  "= (" << beta_params.alpha << ", "<< beta_params.beta << ")\n"
          << std::setw(15) << "  (w0, pi0)" 
          <<  "= (" << w0 << ", "<< pi0 << ")\n";
  }
  
  void set_beta_params(const double alpha_eta, const double beta_eta) {
    beta_params.alpha = alpha_eta;
    beta_params.beta = beta_eta;    
  }
  
  
  
  void set_xi_to_random_labels() {
    for (int j=0; j < J; j++) {
      xi[j] = sample_int_vec(L, n[j]);
      // xi[j] = sample_int_vec(1, n[j]);
    }
  }
  
  void set_z_to_random_labels() {
    z = sample_int_vec(K, J);
    // z = sample_int_vec(1, J);
  }
  
  arma::uvec get_xi_freq_over_z(const int k) {
    arma::uvec count1(L, arma::fill::zeros);
    for (int j = 0; j < xi.size(); j++) {
      if (z(j) == k) {
        count1 += get_freq(xi[j], L);
      }    
    }      
    return count1;      
  }
  
  void update_w(){
    for (int k = 0; k < K; k++) {
      xi_freq_over_z.col(k) = get_xi_freq_over_z(k);
      v.col(k) = gem_gibbs_update_v3(xi_freq_over_z.col(k), w0);
      w.col(k) = stick_break( v.col(k) );
    }        
  }
  
  void update_pi() {
    // pi = stick_break( gem_gibbs_update_v2(get_freq(z, K), pi0) );
    z_freq = get_freq(z, K);
    u = gem_gibbs_update_v3(z_freq, pi0);
    pi = stick_break( u );
  }
  
  void update_z_element(const int j) {
    int zj_old = z(j);
    
    // calculate log_prob for updating z(j)
    List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
    arma::mat D = out["lambda"]; 
    arma::mat Dbar = Rcpp::as<arma::umat>(out["NN"]) - D;
    
    arma::uvec xi_j_freq = get_freq(xi[j], L);
    
    // Rcpp::print(wrap(m.slice(zj_old)-D));
    // Rcpp::print(wrap(mbar.slice(zj_old)-Dbar));
    
    arma::vec log_prob = comp_tensor_log_beta_ratio_sums(
      m, mbar, D, Dbar, zj_old, beta_params.alpha, beta_params.beta
    );
    log_prob += log(w.t() + perturb) * xi_j_freq;
    //for (int ll = 0; ll < L; ll++) { log_prob += log(w[ll, z(j)] + perturb) * xi_j_freq[ll]; }
    
    log_prob += log(pi + perturb);
    
    // update z(j)
    z(j) = sample_index(safe_exp(log_prob)); 
    
    // update m and mbar tensors
    if (z(j) != zj_old) {
      m.slice(z(j)) += D;
      m.slice(zj_old) -= D;
      mbar.slice(z(j)) += Dbar;
      mbar.slice(zj_old) -= Dbar;
      z_freq(z(j))++;
      z_freq(zj_old)--;
    }
  }
  
  void update_xi_element(const int j, const int s) {
    
    arma::vec U = sp_single_col_compress(A[j], s, xi[j], L);
    //  arma::vec U  = blk_compressions[j].row(s).t();
    // Rcpp::print(wrap(U));
    // Rcpp::print(wrap(sp_single_col_compress(A[j], s, xi[j], L)));
    
    arma::uvec V = get_freq(xi[j], L);
    V(xi[j](s))--;
    
    int xi_j_s_old = xi[j](s);
    // prob is K x 1 vector`
    arma::vec temp = comp_log_beta_ratio_sums(m.slice(z(j)), mbar.slice(z(j)), U, V, xi_j_s_old, beta_params.alpha, beta_params.beta);
    
    arma::vec log_prob = temp + log(w.col(z(j)));
    
    xi[j](s) = sample_index(safe_exp(log_prob)); // update z(s) -- this the zs_new we pick
    
    // update m and mbar
    if ( xi_j_s_old != xi[j](s) ) {
      arma::mat D = comp_blk_sums_diff_v1(U, xi[j](s), xi_j_s_old);
      arma::mat DN = comp_blk_sums_diff_v1(arma::conv_to<arma::vec>::from(V), xi[j](s), xi_j_s_old);
      m.slice(z(j)) += D;
      mbar.slice(z(j)) += DN - D;
      // update_col_compress(blk_compressions[j], A[j], s, xi_j_s_old, xi[j](s));
    }
    
    
    // sbm_update_labels(A[j], s, xi[j], L, 
    //             m.slice(z(j)), mbar.slice(z(j)), 
    //             w.col(z(j)), beta_params.alpha, beta_params.beta);
  }
  
  void update_xi_element_marginal(const int j, const int s) {
    
    arma::vec U = sp_single_col_compress(A[j], s, xi[j], L);
    
    arma::uvec V = get_freq(xi[j], L);
    V(xi[j](s))--;
    
    int xi_j_s_old = xi[j](s);
    
    List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
    arma::mat lambda = out["lambda"];
    arma::umat NN = out["NN"];
    
    arma::vec prob(L, arma::fill::zeros);
    for (int k = 0; k < K; k++) {
      arma::vec temp = comp_beta_ratio_prods_v1(m.slice(k) + lambda, mbar.slice(k) + NN - lambda, U, V, xi_j_s_old, beta_params.alpha, beta_params.beta);
      for (int l = 0; l < L; l++) {
        prob(l) += w.col(k)(l) * temp(l) * pi(k);
      }
    }
    
    xi[j](s) = sample_index(prob);
    
    // update m and mbar
    if ( xi_j_s_old != xi[j](s) ) {
      arma::mat D = comp_blk_sums_diff_v1(U, xi[j](s), xi_j_s_old);
      arma::mat DN = comp_blk_sums_diff_v1(arma::conv_to<arma::vec>::from(V), xi[j](s), xi_j_s_old);
      m.slice(z(j)) += D;
      mbar.slice(z(j)) += DN - D;
    }
    
  }
  
  List run_gibbs(const int niter, const int version) {
    // Run full Gibbs updates for "niter" iterations and record label history
    std::vector<std::vector<arma::uvec>> xi_hist(niter+1);
    arma::umat z_hist(J, niter+1);
    
    // xi_hist[0] = xi;
    z_hist.col(0) = z + 1;
    // if (init_count_tensors)
    comp_count_tensors();
    
    set_xi_to_dpsbm_labels(50); // causes an error when updating z
    xi_hist[0] = xi;
    
    for (int iter = 0; iter < niter; iter++) {
      
      // z_freq = get_freq(z, K);
      arma::uvec idx = arma::find(z_freq > 0);
      for (auto k : idx) {
        arma::uvec z_cluster_k = arma::find(z == k);
        for (auto j : z_cluster_k) {
          for (int s = 0; s < n(j); s++) {
            
            switch (version) {
            case 1:
              update_xi_element(j, s);
              break;
            case 2:
              update_xi_element_marginal(j, s);
              break;
            }
          } // s                        
          update_z_element(j);
        }
      }
      
      // for (int j = 0; j < J; j++) {
      //     blk_compressions[j] = sp_compress_col(A[j], xi[j], L);
      // }
      
      update_pi();
      update_w();
      
      // Rcpp::print(wrap(z.t()));
      xi_hist[iter+1] = xi;
      z_hist.col(iter+1) = z + 1;
    } // iter
    
    return Rcpp::List::create( 
      Rcpp::Named("z") = z_hist,
      Rcpp::Named("xi") = xi_hist
    );
  }
  
  void update_pi0(){
    
  }
  
  void update_w0() {
    
  }
  
  // Naive Gibbs ------->
  
  void update_z_element_naive(const int j) {
    // arma::cube  m_old = m;
    // arma::cube  mbar_old = mbar;
    arma::vec log_prob(K, arma::fill::zeros);
    
    // int zj_old = z(j);
    arma::uvec xi_j_freq = get_freq(xi[j], L);
    
    for (int r = 0; r < K; r++) {
      z(j) = r;
      comp_count_tensors(); // this updates "m" and "mbar" based on the current "z"
      arma::cube temp = 
        cube_lbeta(m + beta_params.alpha, mbar + beta_params.beta); // - cube_lbeta(m_old + beta_params.alpha, mbar_old + beta_params.beta);
      
      for (int k = 0; k < K; k++){
        log_prob(r) += arma::sum( arma::trimatu(temp.slice(k)).as_col() );
      }
    }
    
    log_prob += log(w.t() + perturb) * xi_j_freq;
    log_prob += log(pi + perturb);
    
    // update z(j)
    z(j) = sample_index(safe_exp(log_prob)); 
    comp_count_tensors();
  }
  
  List run_gibbs_naive(const int niter, const int version) {
    // Run full Gibbs updates for "niter" iterations and record label history
    std::vector<std::vector<arma::uvec>> xi_hist(niter+1);
    arma::umat z_hist(J, niter+1);
    
    set_xi_to_dpsbm_labels(50);
    xi_hist[0] = xi;
    z_hist.col(0) = z + 1;
    // if (init_count_tensors) 
    comp_count_tensors();
        
    for (int iter = 0; iter < niter; iter++) {
      for (int j = 0; j < J; j++) {
        update_z_element_naive(j);
        for (int s = 0; s < n(j); s++) {
          switch (version) {
          case 1:
            update_xi_element(j, s);
            break;
          case 2:
            update_xi_element_marginal(j, s);
            break;
          }
        } // s
        //update_z_element_naive(j); // error if we run this
      } // j
      update_w();
      update_pi();
      xi_hist[iter+1] = xi;
      z_hist.col(iter+1) = z + 1;
    } // iter
    
    return Rcpp::List::create( 
      Rcpp::Named("z") = z_hist,
      Rcpp::Named("xi") = xi_hist
    );
  }
  
  // end Naive Gibbs ------->
  
  // void update_xi_element_naive(const int j, const int s) {
  //     arma::vec log_prob(L, arma::fill::zeros);
  
  //     for (int ll = 0; ll < L; ll++) {
  //         xi[j](s) = ll;
  //         comp_count_tensors(); // this updates "m" and "mbar" based on the current "z"
  //         arma::cube temp = 
  //             cube_lbeta(m + beta_params.alpha, mbar + beta_params.beta); // - cube_lbeta(m_old + beta_params.alpha, mbar_old + beta_params.beta);
  
  //         for (int k = 0; k < K; k++){
  //             log_prob(ll) += arma::sum( arma::trimatu(temp.slice(k)).as_col() );
  //         }
  //     }
  
  //     log_prob += log(w.col(z(j)) + perturb);
  
  //     // update xi[j](s)
  //     xi[j](s) = sample_index(safe_exp(log_prob)); 
  //     comp_count_tensors();
  // }
  
private:
  const double perturb = 1e-11;
};


// void test_nsbm_cpp(List A, const int K, const int L) {

//     NestedSBM mynsbm(A, K, L);

//     mynsbm.run_gibbs(100);
// }

RCPP_MODULE(sbm_module) {
  class_<NestedSBM>("NestedSBM")
  .constructor<List, int, int>()
  .field("A", &NestedSBM::A)
  .field("J", &NestedSBM::J)
  .field("K", &NestedSBM::K)
  .field("L", &NestedSBM::L)
  .field("n", &NestedSBM::n)
  .field("z", &NestedSBM::z)
  .field("xi", &NestedSBM::xi)
  .field("m", &NestedSBM::m)
  .field("mbar", &NestedSBM::mbar)
  .field("w", &NestedSBM::w)
  .field("pi", &NestedSBM::pi)
  .field("w0", &NestedSBM::w0)
  .field("pi0", &NestedSBM::pi0)
  .field("xi_freq_over_z", &NestedSBM::xi_freq_over_z)
  .field("z_freq", &NestedSBM::z_freq)
  .method("set_beta_params", &NestedSBM::set_beta_params)
  .method("comp_count_tensors", &NestedSBM::comp_count_tensors)
  .method("print", &NestedSBM::print)
  .method("update_z_element", &NestedSBM::update_z_element)
  .method("update_z_element_naive", &NestedSBM::update_z_element_naive)
  .method("update_xi_element", &NestedSBM::update_xi_element)
  .method("set_xi_to_random_labels", &NestedSBM::set_xi_to_random_labels)
  .method("set_z_to_random_labels", &NestedSBM::set_z_to_random_labels)
  .method("get_xi_freq_over_z", &NestedSBM::get_xi_freq_over_z)
  .method("run_gibbs", &NestedSBM::run_gibbs)
  .method("run_gibbs_naive", &NestedSBM::run_gibbs_naive)
  .method("run_gibbs_via_eta", &NestedSBM::run_gibbs_via_eta)
  .method("update_w", &NestedSBM::update_w)
  .method("update_pi", &NestedSBM::update_pi)
  ;
};

RCPP_EXPOSED_CLASS(NestedSBM);