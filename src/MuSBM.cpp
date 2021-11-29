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


class MultSBM {
    public:
        // List A;
        std::vector<arma::sp_mat> A;
        int J;  // number of networks
        int L;  // truncation level for the prior on "w"
        arma::uvec n; // n(j) = number nodes in network j
        std::vector<arma::uvec> xi; // xi[j](s) = label of nodes s in network j
        arma::vec w;    // the prior on "xi"

        arma::mat u;
        arma::mat v;
        arma::mat eta;
        
        BetaParameters beta_params;
        // double w0;
        // double pi0;

        // arma::vec z_log_prob_record; // for diagnostics

        MultSBM(// const List A_,
            const std::vector<arma::sp_mat>& A_, 
            // const arma::uvec z_init, 
            const int L,
            const double alpha_eta, 
            const double beta_eta) : L{L} {

            // initialize and allocate variables
            A = A_;
            //J = A.length();
            J = A.size();
            n = arma::uvec(J);
            xi = std::vector<arma::uvec>(J);

             for (int j=0; j < J; j++) {
                //n(j) = Rcpp::as<arma::sp_mat>(A[j]).n_rows;
                n(j) = A[j].n_rows;
            }

            beta_params.alpha = alpha_eta;
            beta_params.beta = beta_eta;
            set_xi_to_random_labels();

            eta = arma::mat(L, L, arma::fill::zeros);   
            u = arma::mat(L, L, arma::fill::zeros);
            v = arma::mat(L, L, arma::fill::zeros);
            w = arma::vec(L, arma::fill::zeros);

            // Rcpp::print(wrap( blk_compressions[1]));
        }

        void set_beta_params(const double alpha_eta, const double beta_eta) {
            beta_params.alpha = alpha_eta;
            beta_params.beta = beta_eta;    
        }

        void set_xi_to_random_labels() {
            for (int j=0; j < J; j++) {
                xi[j] = sample_int_vec(L, n[j]);
            }
        }

        void update_eta() {    
            // Update the eta-related tensors
            arma::mat m(L, L, arma::fill::zeros); 
            arma::mat mbar(L, L, arma::fill::zeros);

            for (int j=0; j < J; j++) {
                List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
                arma::mat lambda = out["lambda"];
                arma::umat NN = out["NN"]; 
                m += lambda;
                mbar += NN - lambda;
            }

            eta = symmat_rbeta(m + beta_params.alpha, mbar + beta_params.beta);

            u = log(eta/(1-eta) + perturb);
            v = log(1-eta + perturb);
        }

        void update_w() {
            arma::vec nn(L, arma::fill::zeros);
            for (int j = 0; j < xi.size(); j++) {
                nn += arma::conv_to<arma::vec>::from(
                    get_freq(xi[j], L)
                );
            }     
            w = rdirichlet(nn + 1);
        }

        void update_xi_element(const int j, const int s) {
            
            arma::vec taui = sp_single_col_compress(A[j], s, xi[j], L);
            arma::uvec mmi = get_freq(xi[j], L);
            mmi(xi[j](s))--;

            arma::vec log_prob = u * taui + v * mmi + log(w);

            xi[j](s) = sample_index(safe_exp(log_prob));
        }
       
       std::vector<std::vector<arma::uvec>> run_gibbs(const int niter) {
            // Run full Gibbs updates for "niter" iterations and record label history
            
           std::vector<std::vector<arma::uvec>> xi_hist(niter+1);
           xi_hist[0] = xi;       
            
            // comp_count_tensors();
            for (int iter = 0; iter < niter; iter++) {   
                update_eta();
                update_w();
                for (int j = 0; j < J; j++) {
                    for (int s = 0; s < n(j); s++) {
                        update_xi_element(j, s);
                    } // s
                } // j
                // Rcpp::print(wrap(z.t()));
                xi_hist[iter+1] = xi;
            } // iter

            return xi_hist;
            // return Rcpp::List::create( 
            //     Rcpp::Named("z") = z_hist,
            //     Rcpp::Named("xi") = xi_hist
            // );
        }

    private:
        const double perturb = 1e-11;
};


// void test_nsbm_cpp(List A, const int K, const int L) {

//     MultSBM mynsbm(A, K, L);

//     mynsbm.run_gibbs(100);
// }

RCPP_MODULE(sbm_module) {
      class_<MultSBM>("MultSBM")
      // .constructor<List, int, double, double>()
      .constructor<std::vector<arma::sp_mat>, int, double, double>()
      .field("A", &MultSBM::A)
      .field("J", &MultSBM::J)
      .field("L", &MultSBM::L)
      .field("n", &MultSBM::n)
      .field("xi", &MultSBM::xi)
      .field("eta", &MultSBM::eta)
      .method("set_beta_params", &MultSBM::set_beta_params)
//      .method("print", &MultSBM::print)
      .method("update_xi_element", &MultSBM::update_xi_element)
      .method("set_xi_to_random_labels", &MultSBM::set_xi_to_random_labels)
      .method("run_gibbs", &MultSBM::run_gibbs)
      .method("update_eta", &MultSBM::update_eta)
      .method("update_w", &MultSBM::update_w)
      ;
};

RCPP_EXPOSED_CLASS(MultSBM);