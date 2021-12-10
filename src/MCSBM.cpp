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

    BetaParameters& operator=(const BetaParameters& that)  {
        alpha = that.alpha;
        beta = that.beta;
        return *this;
    }

};


class MCSBM {
    public:
        // List A;
        std::vector<arma::sp_mat> A;
        int J;  // number of networks
        int K;  // truncation level for the prior on "z"
        int L;  // truncation level for the prior on "w"
        arma::uvec n; // n(j) = number nodes in network j
        arma::uvec z; // z(j) = label of network j
        std::vector<arma::uvec> xi; // xi[j](s) = label of nodes s in network j
        
        
        arma::mat w;    // the prior on "xi"
        arma::vec pi;   // the prior on "z"

        arma::cube u;
        arma::cube v;
        arma::cube eta;

        arma::cube m;   // m(x,y,j)
        arma::ucube N; // N(x,y,j)
        
        BetaParameters beta_params;
        double rnd_prob;  // Random z-label probability
        double sa_temp;  // Simulated Annealing Temperature 
        double decay;
        // double w0;
        // double pi0;

        // arma::vec z_log_prob_record; // for diagnostics

        MCSBM(const std::vector<arma::sp_mat>& A_,
            const int K,
            const int L,
            const double alpha_eta, 
            const double beta_eta) : L{L}, K{K} { //, rnd_prob{rnd_prob}, decay{decay} {

            // initialize and allocate variables
            A = A_;
            //J = A.length();
            J = A.size();
            n = arma::uvec(J);
            xi = std::vector<arma::uvec>(J);
            sa_temp = 1;
            rnd_prob = 0;
            decay = 1;

            for (int j=0; j < J; j++) {
                //n(j) = Rcpp::as<arma::sp_mat>(A[j]).n_rows;
                n(j) = A[j].n_rows;
            }

            beta_params.alpha = alpha_eta;
            beta_params.beta = beta_eta;
            set_xi_to_random_labels();
            set_z_to_random_labels();

            eta = arma::cube(L, L, K, arma::fill::zeros);   
            u = arma::cube(L, L, K, arma::fill::zeros);
            v = arma::cube(L, L, K, arma::fill::zeros);
            w = arma::mat(L, K, arma::fill::zeros);

            m = arma::cube(L, L, J, arma::fill::zeros);  // m tensor
            N = arma::ucube(L, L, J, arma::fill::zeros);  // N tensor

            // Rcpp::print(wrap( blk_compressions[1]));
        }

        // copy assignment operator 
        // MCSBM& operator=(const MCSBM&) = default;

        void set_beta_params(const double alpha_eta, const double beta_eta) {
            beta_params.alpha = alpha_eta;
            beta_params.beta = beta_eta;    
        }

        void set_xi_to_given(std::vector<arma::uvec> xi_given) {
            for (int j=0; j < J; j++) {
                xi[j] = xi_given[j]-1;
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

        void update_count_tensors() {
            for (int j=0; j < J; j++) {
                update_count_tensor(j);
            }
        }

        void update_count_tensor(const int j) {
            List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
            arma::mat lambda = out["lambda"];
            arma::umat NN = out["NN"]; 
            m.slice(j) += lambda;
            N.slice(j) += NN;
        }

        // Assumes that m and N are up-to-date
        void update_eta() {  
            // Update the eta-related tensors: eta, u and v
            // comp_count_tensors();
            arma::cube m_sum(L, L, K, arma::fill::zeros);
            arma::ucube N_sum(L, L, K, arma::fill::zeros);
            for (int j=0; j < J; j++) {
                m_sum.slice(z(j)) += m.slice(j);
                N_sum.slice(z(j)) += N.slice(j);
            }

            for (int k = 0; k < K; k++) {
                eta.slice(k) = symmat_rbeta(m_sum.slice(k) + beta_params.alpha, 
                                            N_sum.slice(k) - m_sum.slice(k) + beta_params.beta);
            }  
            u = log(eta/(1-eta) + perturb);
            v = log(1-eta + perturb);
        }

        void update_w() {
            for (int k = 0; k < K; k++) {
                arma::vec nn(L, arma::fill::zeros);
                for (int j = 0; j < xi.size(); j++) {
                    if (z(j) == k) {
                        nn += arma::conv_to<arma::vec>::from( get_freq(xi[j], L) );
                    }  
                }  
                w.col(k) = rdirichlet(nn + 1);
            }  
            
            // arma::vec nn(L, K, arma::fill::zeros);
            // for (int j = 0; j < xi.size(); j++) {
            //     nn.col(z(j)) += arma::conv_to<arma::vec>::from(
            //         get_freq(xi[j], L)
            //     );
            // }
            // Rcout << "1\n";
            // for (int k = 0; k < K; k++) {
            //     w.col(k) = rdirichlet(nn.col(k) + 1);
            // }
        }

         void update_pi() {
            arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
            pi = rdirichlet(nn + 1);
        }
         
         // Assumes update m and N tensors
         void update_z_element(const int j, const int alt_value) {

            if (R::runif(0,1) > rnd_prob) { // update according to Gibbs dynamic
                arma::vec log_prob(K, arma::fill::zeros);

                // int zj_old = z(j);
                arma::uvec xi_j_freq = get_freq(xi[j], L);

                for (int r = 0; r < K; r++) {
                    for (int x = 0; x < L; x++) {
                        for (int y = x; y < L; y++) {
                            log_prob(r) += u(x,y,r) * m(x,y,j) + v(x,y,r)* N(x,y,j);
                        }
                    }
                    // log_prob(r) += arma::sum( arma::trimatu(temp.slice(k)).as_col() );
                }
                
                log_prob += log(w.t() + perturb) * xi_j_freq;
                log_prob += log(pi + perturb);

                // update z(j)
            
                z(j) = sample_index(safe_exp(log_prob / sa_temp)); 
            } else { // update using the provided alt_value
                z[j] = alt_value;
            }
            // comp_count_tensors();
        }


        void update_xi_element(const int j, const int s) {
            
            arma::vec taui = sp_single_col_compress(A[j], s, xi[j], L);
            arma::uvec mmi = get_freq(xi[j], L);
            mmi(xi[j](s))--;

            arma::vec log_prob = u.slice(z(j)) * taui + v.slice(z(j)) * mmi + log(w.col(z(j)));

            xi[j](s) = sample_index(safe_exp(log_prob));
        }

        void run_gibbs_step(const arma::uvec& alt_vec) {

            update_count_tensors();
            update_eta();
            update_w();
            update_pi();
            for (int j = 0; j < J; j++) {
                for (int s = 0; s < n(j); s++) {
                    update_xi_element(j, s);
                } // s             
                update_z_element(j, alt_vec[j]);
            } // j
        }

        void update_rnd_prob() {
             rnd_prob *= decay;  
        }

        void update_sa_temp() {
            sa_temp = std::max(sa_temp*decay, 1.0);
        }

       // std::vector<std::vector<arma::uvec>> 
       List run_gibbs(const int niter) {
            // Run full Gibbs updates for "niter" iterations and record label history
            
            std::vector<std::vector<arma::uvec>> xi_hist(niter+1);
            arma::umat z_hist(J, niter+1);
            xi_hist[0] = xi;
            z_hist.col(0) = z + 1;       
            
            // update_count_tensors(); 
            for (int iter = 0; iter < niter; iter++) {
                run_gibbs_step(sample_int_vec(K, J)); 
                update_sa_temp();
                update_rnd_prob();
              
                // Rcpp::print(wrap(z.t()));
                z_hist.col(iter+1) = z + 1;
                xi_hist[iter+1] = xi;
            } // iter

            // return xi_hist;
            return Rcpp::List::create( 
                 Rcpp::Named("z") = z_hist,
                 Rcpp::Named("xi") = xi_hist
            );
        }

    private:
        const double perturb = 1e-11;
};

// [[Rcpp::export]]
List mix_mcsbm(const std::vector<arma::sp_mat>& A, 
            const int K,
            const int L,
            const double alpha_eta, 
            const double beta_eta,
            const double decay, 
            const double sa_temp, 
            const double rnd_prob, 
            const int niter,
            const int n_models = 3) {
    
    
    std::vector<MCSBM> models(n_models, MCSBM(A, K, L, alpha_eta, beta_eta));
  
    for (auto &model: models) {
        model.rnd_prob = rnd_prob;
        model.decay = decay;
        model.sa_temp = sa_temp;
        model.set_xi_to_random_labels();
        model.set_z_to_random_labels();
    }
  
    // MCSBM mod1(A, K, L, alpha_eta, beta_eta, rnd_prob, decay);
    // MCSBM mod2(A, K, L, alpha_eta, beta_eta, rnd_prob, decay);
    // xi_hist[0] = mod1.xi;
    // z_hist.col(0) = mod1.z + 1;   


    int J = A.size();
    std::vector<std::vector<arma::uvec>> xi_hist(niter+1);
    arma::umat z_hist(J, niter+1);

    xi_hist[0] = models[0].xi;
    z_hist.col(0) = models[0].z + 1;       
    

    auto sample_except_r = [](const int n_models, const int r) {
        int x = sample_int(n_models-1);
        if (x >= r) x++;
        return x;
    };
    

    for (int iter = 0; iter < niter; iter++) {
        // mod1.run_gibbs_step(mod2.z);
        // mod2.run_gibbs_step(mod1.z);
        // models[0].run_gibbs_step(models[1].z);
        // models[1].run_gibbs_step(models[0].z);
        for (int r = 0; r < n_models; r++) {
        //    arma::uvec alt_z(J);
        //    for (int j =0; j < J; j++) {
        //        int alt_model_index = sample_except_r(n_models, r);
        //        alt_z[j] = models[alt_model_index].z[j];
        //    } 
        //    models[r].run_gibbs_step(alt_z);

           int alt_model_index = sample_except_r(n_models, r); 
           models[r].run_gibbs_step(models[alt_model_index].z);
           models[r].update_sa_temp();
           models[r].update_rnd_prob();
           // models[r].rnd_prob *= decay;
        }        
        
        z_hist.col(iter+1) = models[0].z + 1;
        xi_hist[iter+1] = models[0].xi;
        // z_hist.col(iter+1) = mod1.z + 1;
        // xi_hist[iter+1] = mod1.xi;
    } // iter

    // return xi_hist;
    return Rcpp::List::create( 
            Rcpp::Named("z") = z_hist,
            Rcpp::Named("xi") = xi_hist
    );
    
}

RCPP_MODULE(sbm_module) {
      class_<MCSBM>("MCSBM")
      .constructor<std::vector<arma::sp_mat>, int, int, double, double>()
      .field("A", &MCSBM::A)
      .field("J", &MCSBM::J)
      .field("K", &MCSBM::K)
      .field("L", &MCSBM::L)
      .field("n", &MCSBM::n)
      .field("xi", &MCSBM::xi)
      .field("z", &MCSBM::z)
      .field("w", &MCSBM::w)
      .field("pi", &MCSBM::pi)
      .field("eta", &MCSBM::eta)
      .field("m", &MCSBM::m)
      .field("N", &MCSBM::N)
      .field("rnd_prob", &MCSBM::rnd_prob)
      .field("sa_temp", &MCSBM::sa_temp)
      .field("decay", &MCSBM::decay)
      .method("set_beta_params", &MCSBM::set_beta_params)
//      .method("print", &MCSBM::print)
      .method("update_xi_element", &MCSBM::update_xi_element)
      .method("update_z_element", &MCSBM::update_z_element)
      .method("set_z_to_random_labels", &MCSBM::set_z_to_random_labels)
      .method("set_xi_to_random_labels", &MCSBM::set_xi_to_random_labels)
      .method("run_gibbs", &MCSBM::run_gibbs)
      .method("update_eta", &MCSBM::update_eta)
      .method("update_w", &MCSBM::update_w)
      .method("update_pi", &MCSBM::update_pi)
      .method("set_xi_to_given", &MCSBM::set_xi_to_given)
      .method("update_count_tensors", &MCSBM::update_count_tensors)
      ;
};

RCPP_EXPOSED_CLASS(MCSBM);