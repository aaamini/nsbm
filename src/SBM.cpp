// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "utils.h"
#include "sampling.h"
// #include "beta_calcs.h"
// #include "dpsbm.h"

using namespace Rcpp;

struct BetaParameters {
    double alpha;
    double beta;

    BetaParameters() :alpha{1}, beta{1} {};
    BetaParameters(const double a, const double b) :alpha{a}, beta{b} {};
};


class SBM {
    public:
        arma::sp_mat A;
        int K;
        int n;
        arma::uvec z;
        arma::mat eta;
        arma::mat u;
        arma::mat v;
        arma::vec pri; 
        
        BetaParameters beta_params;
        // double w0;
        // double pi0;

        // arma::vec z_log_prob_record; // for diagnostics

        SBM(const arma::sp_mat A_, 
            // const arma::uvec z_init, 
            const int K,
            const double alpha_eta, 
            const double beta_eta) : K{K} {

            // initialize and allocate variables
            A = A_;
            n = A.n_rows;
            beta_params.alpha = alpha_eta;
            beta_params.beta = beta_eta;
            set_z_to_random_labels();

            eta = arma::mat(K, K, arma::fill::zeros);   
            u = arma::mat(K, K, arma::fill::zeros);
            v = arma::mat(K, K, arma::fill::zeros);
            pri = arma::vec(K, arma::fill::zeros);

            // Rcpp::print(wrap( blk_compressions[1]));
        }

        void set_beta_params(const double alpha_eta, const double beta_eta) {
            beta_params.alpha = alpha_eta;
            beta_params.beta = beta_eta;    
        }

        void set_z_to_random_labels() {
             z = sample_int_vec(K, n);
        }

        void update_eta() {    
            // Update the eta-related tensors
            List out = comp_blk_sums_and_sizes(A, z, K);
            arma::mat lambda = out["lambda"];
            arma::umat NN = out["NN"]; 
            eta = symmat_rbeta(lambda + beta_params.alpha, NN - lambda + beta_params.beta);

            u = log(eta/(1-eta) + perturb);
            v = log(1-eta + perturb);
        }

        void update_pri() {
            
            arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
            pri = rdirichlet(nn + 1);
        }

        void update_z_element(const int i) {

            arma::vec taui = sp_single_col_compress(A, i, z, K);
            arma::uvec mmi = get_freq(z, K);
            mmi(z(i))--;
            arma::vec log_prob = u * taui + v * mmi + log(pri);

            z(i) = sample_index(safe_exp(log_prob));
        }
        
       
        arma::umat run_gibbs(const int niter) {
            // Run full Gibbs updates for "niter" iterations and record label history
            
            arma::umat z_hist(n, niter+1);
            z_hist.col(0) = z + 1;       
            
            // comp_count_tensors();
            for (int iter = 0; iter < niter; iter++) {   
                update_eta();
                update_pri();
                for (int i = 0; i < n; i++) {
                    update_z_element(i);
                }
                // Rcpp::print(wrap(z.t()));
                z_hist.col(iter+1) = z + 1;
            } // iter

            return z_hist;
            // return Rcpp::List::create( 
            //     Rcpp::Named("z") = z_hist,
            //     Rcpp::Named("xi") = xi_hist
            // );
        }

    private:
        const double perturb = 1e-11;
};


// void test_nsbm_cpp(List A, const int K, const int L) {

//     SBM mynsbm(A, K, L);

//     mynsbm.run_gibbs(100);
// }

RCPP_MODULE(sbm_module) {
      class_<SBM>("SBM")
      .constructor<arma::sp_mat, int, double, double>()
      .field("A", &SBM::A)
      .field("K", &SBM::K)
      .field("n", &SBM::n)
      .field("z", &SBM::z)
      .field("eta", &SBM::eta)
      .method("set_beta_params", &SBM::set_beta_params)
//      .method("print", &SBM::print)
      .method("update_z_element", &SBM::update_z_element)
      .method("set_z_to_random_labels", &SBM::set_z_to_random_labels)
      .method("run_gibbs", &SBM::run_gibbs)
      .method("update_eta", &SBM::update_eta)
      .method("update_pri", &SBM::update_pri)
      ;
};

RCPP_EXPOSED_CLASS(SBM);