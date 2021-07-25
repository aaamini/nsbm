// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "utils.h"
#include "sampling.h"
#include "beta_calcs.h"

using namespace Rcpp;

struct BetaParameters {
    double alpha;
    double beta;

    BetaParameters() :alpha{1}, beta{1} {};
    BetaParameters(const double a, const double b) :alpha{a}, beta{b} {};
};

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
        arma::cube m;   // m(x,y,k)
        arma::cube mbar; // mbar(x,y,k)
        BetaParameters beta_params;

        NestedSBM(List A_, 
            // const arma::uvec z_init, 
            const int K, const int L, 
            const double alpha_eta = 1, 
            const double beta_eta = 1,
            const double w0 = 1, 
            const double pi0 = 1) : K{K}, L{L}, w0{w0}, pi0{pi0} {

            // initialize and allocate variables
            A = A_;
            J = A.length();
            n = arma::uvec(J);
            xi = std::vector<arma::uvec>(J);
            // z = arma::uvec(J);
            m = arma::cube(L, L, K, arma::fill::zeros); // m tensor
            mbar = arma::cube(L, L, K, arma::fill::zeros); // mbar tensor
            pi = arma::vec(K, arma::fill::ones);
            w = arma::mat(L, K, arma::fill::ones);

            for (int j=0; j < J; j++) {
                n(j) = Rcpp::as<arma::sp_mat>(A[j]).n_rows;
                // xi[j] = sample_int_vec(L, n[j]);
            }
            // z = sample_int_vec(K, J);
            beta_params.alpha = alpha_eta;
            beta_params.beta = beta_eta;
            set_xi_to_random_labels();
            set_z_to_random_labels();
        }

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

        void set_xi_to_random_labels() {
            for (int j=0; j < J; j++) {
                xi[j] = sample_int_vec(L, n[j]);
            }
        }

        void set_z_to_random_labels() {
            z = sample_int_vec(K, J);
        }

        void update_z_element(const int j) {
            int zj_old = z(j);

            // calculate log_prob for updating z(j)
            List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
            arma::mat D = out["lambda"]; 
            arma::mat Dbar = Rcpp::as<arma::umat>(out["NN"]) - D;

            arma::uvec xi_j_freq = get_freq(xi[j], L);
            arma::vec log_prob = comp_tensor_log_beta_ratio_sums(
                m, mbar, D, Dbar, zj_old, beta_params.alpha, beta_params.beta
            );
            log_prob += log(w.t() + perturb) * xi_j_freq;
            log_prob += log(pi + perturb);

            // update z(j)
            z(j) = sample_index(safe_exp(log_prob)); 

            // update m and mbar tensors
            if (z(j) != zj_old) {
                m.slice(z(j)) += D;
                m.slice(zj_old) -= D;
                mbar.slice(z(j)) += Dbar;
                mbar.slice(zj_old) -= Dbar;
            }
        }

        void update_xi_element(const int j, const int s) {

            sbm_update_labels(A[j], s, xi[j], L, 
                        m.slice(z(j)), mbar.slice(z(j)), 
                        w.col(z(j)), beta_params.alpha, beta_params.beta);
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
                w.col(k) = stick_break( gem_gibbs_update_v2(get_xi_freq_over_z(k), w0) );
            }        
        }

        void update_pi() {
            pi = stick_break( gem_gibbs_update_v2(get_freq(z, K), pi0) );
        }

        // List run_gibbs(const int niter = 100, const bool init_count_tensors = true) {

        List run_gibbs(const int niter) {
            // Run full Gibbs updates for "niter" iterations and record label history
            std::vector<std::vector<arma::uvec>> xi_hist(niter);
            arma::umat z_hist(J, niter);

            // if (init_count_tensors) 
            comp_count_tensors();

            for (int iter = 0; iter < niter; iter++) {
                update_w();
                update_pi();
                for (int j = 0; j < J; j++) {
                    for (int s = 0; s < n(j); s++) {
                        update_xi_element(j, s);
                    } // s
                    update_z_element(j);
                } // j
                xi_hist[iter] = xi;
                z_hist.col(iter) = z + 1;
            } // iter

            return Rcpp::List::create( 
                Rcpp::Named("z") = z_hist,
                Rcpp::Named("xi") = xi_hist
            );
        }

    private:
        const double perturb = 1e-11;
        const double w0;
        const double pi0;
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
      .method("comp_count_tensors", &NestedSBM::comp_count_tensors)
      .method("print", &NestedSBM::print)
      .method("update_z_element", &NestedSBM::update_z_element)
      .method("update_xi_element", &NestedSBM::update_xi_element)
      .method("set_xi_to_random_labels", &NestedSBM::set_xi_to_random_labels)
      .method("set_z_to_random_labels", &NestedSBM::set_z_to_random_labels)
      .method("get_xi_freq_over_z", &NestedSBM::get_xi_freq_over_z)
      .method("run_gibbs", &NestedSBM::run_gibbs)
      .method("update_w", &NestedSBM::update_w)
      .method("update_pi", &NestedSBM::update_pi)
      ;
};

RCPP_EXPOSED_CLASS(NestedSBM);