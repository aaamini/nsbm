// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

#include "sampling.h" 
#include "utils.h" 

using namespace Rcpp;

// [[Rcpp::export]]
arma::cube comp_beta_tensor(List A, std::vector<arma::uvec> & xi, 
                arma::uvec& z, const int L, const int K,
                const double alpha = 1, const double beta = 1) {
    
    arma::cube m_tensor(L, L, K, arma::fill::zeros); // m tensor
    arma::cube mbar_tensor(L, L, K, arma::fill::zeros); // m tensor
    arma::cube beta_tensor(L, L, K, arma::fill::ones); // beta tensor
    
    for (int j=0; j < A.length(); j++) {
        List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
        arma::mat lambda = out["lambda"];
        arma::umat NN = out["NN"]; 
        m_tensor.slice(z(j)) += lambda;
        mbar_tensor.slice(z(j)) += NN - lambda;
    }

    for (int k=0; k < K; k++) {
        beta_tensor.slice(k) = 
            beta_fun_symmat(m_tensor.slice(k) + alpha, 
                            mbar_tensor.slice(k) + beta);
    }
    
    return beta_tensor;
    // return Rcpp::List::create( 
    //     Rcpp::Named("m") = m_tensor,
    //     Rcpp::Named("mbar") = mbar_tensor,
    //     Rcpp::Named("beta") = beta_tensor
    // );
}

// [[Rcpp::export]]
double cube_all_elem_prod(arma::cube tensor) {
    double out = 1;
    for (int k = 0; k < tensor.n_slices; k++) {
        out *= arma::prod(tensor.slice(k).as_col());
    }
    return out;
}

// [[Rcpp::export]]
arma::vec ndp_gem_gibbs_update(std::vector<arma::uvec> & xi, 
                arma::uvec& z, const int k, 
                const int L, const int K, 
                const double concent_param) {
    
    arma::uvec count1(L, arma::fill::zeros);
    for (int j = 0; j < xi.size(); j++) {
        if (z(j) == k) {
             count1 += get_freq(xi[j], L);
        }    
    }
    arma::uvec count2 = get_up_freq(count1);
    // Rcpp::print(Rcpp::wrap(count1));
    // Rcpp::print(Rcpp::wrap(count2));

    return(
        rbeta_vec(arma::conv_to<arma::vec>::from(count1) + 1, 
                    arma::conv_to<arma::vec>::from(count2) + concent_param)
    );
    
}

// [[Rcpp::export]]
arma::umat fit_nsbm(List A,
                double alpha_eta = 1, double beta_eta = 5,
                double pi0 = 1, double w0 = 1,
                int niter = 50, int K=5, int L = 10, bool verb = true,
                const double perturb = 1e-11,
                const bool rand_init = true, const bool seq_g_update = true) {

    int J = A.length();
//   int tmax = A.length();
    arma::uvec n(J);

    for (int j=0; j < J; j++) {
        n(j) = Rcpp::as<arma::sp_mat>(A[j]).n_rows;
    }

    // return n;

    std::vector<arma::uvec> xi(J); // xi(j)[i] is in 0,1,...,L-1 where i ranges in 0,1,...,n(j)-1
    arma::uvec z(J); // z(k) is in 0,1,...,K-1
    arma::umat z_hist(J, niter, arma::fill::zeros);
    arma::vec  pi(K);
    arma::mat  w(L, K);
    // arma::cube beta_tensor(L, L, K, arma::fill::ones);
    
    // initialize
    for (int j=0; j < J; j++) {
        xi[j] = sample_int_vec(L, n[j]);
    }
    z = sample_int_vec(K, J);

    // initial computation of the m-tensor
    arma::cube m(L, L, K, arma::fill::zeros); // m tensor
    arma::cube mbar(L, L, K, arma::fill::zeros); // mbar tensor
    // arma::cube beta_tensor(L, L, K, arma::fill::ones); // beta tensor
   
    for (int j=0; j < A.length(); j++) {
        List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
        arma::mat lambda = out["lambda"];
        arma::umat NN = out["NN"]; 
        m.slice(z(j)) += lambda;
        mbar.slice(z(j)) += NN - lambda;
    }

    // for (int k=0; k < K; k++) {
    //     beta_tensor.slice(k) = 
    //         beta_fun_symmat(m_tensor.slice(k) + alpha_eta, 
    //                         mbar_tensor.slice(k) + beta_eta);
    // }

    // main loop
    for (int iter = 0; iter < niter; iter++) {

         // updating w
        for (int k = 0; k < K; k++) {
            w.col(k) = stick_break( ndp_gem_gibbs_update(xi, z, k, L, K, w0) );
        }        

        // updating pi
        pi = stick_break( gem_gibbs_update(z, K, pi0) );
        
        // updating xi
        for (int j = 0; j < J; j++) {
            for (int s = 0; s < n(j); s++) {
                // beta_tensor = comp_beta_tensor(A, xi, z, L, K, 
                //                         alpha_eta, beta_eta);
                arma::vec U = sp_single_col_compress(A[j], s, xi[j], L);
                arma::uvec V = get_freq(xi[j], L);
                V(xi[j](s))--;

                int xi_s_j_old = xi[j](s);
                arma::vec log_prob = 
                    comp_log_beta_ratio_sums(m.slice(z(j)), mbar.slice(z(j)), U, V, xi_s_j_old, alpha_eta, beta_eta) +
                    log(w.col(z(j)));

                // arma::vec prob(L, arma::fill::ones);
                // for (int rp = 0; rp < L; rp++) { // rp is the potential new value of z(s)
                //     arma::uvec xi_j_new = xi[j];
                //     xi_j_new(s) = rp;
                    
                //      arma::mat new_bet = comp_beta_matrix(A[j], xi_j_new, L, alpha_eta, beta_eta);
                //     prob(rp) *= arma::prod((new_bet / beta_tensor.slice(z(j))).as_col()) * w(rp, z(j));
                //     // prob(rp) *= static_cast<double>((nn(rp) + 1)) / nn(z(s)); 
                // }
                // // Rcpp::Rcout << "---\n" << prob;
                xi[j](s) = sample_index(safe_exp(log_prob)); ; // update z

                arma::mat D = comp_blk_sums_diff_v1(U, xi[j](s), xi_s_j_old);
                arma::mat DN = comp_blk_sums_diff_v1(arma::conv_to<arma::vec>::from(V), xi[j](s), xi_s_j_old);
                m.slice(z(j)) += D;
                mbar.slice(z(j)) += DN - D;

                // List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
                // arma::mat  lambda = out["lambda"];
                // arma::umat NN = out["NN"]; 
                // m_tensor.slice(z(j)) += lambda;
            }
        }

        // updating z
        for (int j = 0; j < J; j++) {
            // beta_tensor = comp_beta_tensor(A, xi, z, L, K, 
            //                         alpha_eta, beta_eta);

            int zj_old = z(j);

            // calculate log_prob for updating z(j)
            List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
            arma::mat D = out["lambda"]; 
            arma::mat Dbar = Rcpp::as<arma::umat>(out["NN"]) - D;

            arma::uvec xi_j_freq = get_freq(xi[j], L);
            arma::vec log_prob = comp_tensor_log_beta_ratio_sums(
                m, mbar, D, Dbar, zj_old, alpha_eta, beta_eta
            );
            // Rcout << "zj_old = " << zj_old << "log_prob1 = " << log_prob << ", ";
            log_prob += log(w.t() + perturb) * xi_j_freq;
            log_prob += log(pi + perturb);
            // Rcout << "log_prob2 = " << log_prob << "\n";

            // arma::vec prob(L, arma::fill::ones);
            // for (int rp = 0; rp < K; rp++) { // rp is the potential new value of z(s)
            //     arma::uvec z_new = z;
            //     z_new(j) = rp;
            //     arma::cube new_beta_tensor = 
            //         comp_beta_tensor(A, xi, z_new, L, K, 
            //                         alpha_eta, beta_eta);

            //     prob(rp) *= cube_all_elem_prod(new_beta_tensor / beta_tensor);
            //     prob(rp) *= exp( sum(xi_j_freq % log(w.col(rp) + perturb)) );
            //     prob(rp) *= pi(rp);
            //     // prob(rp) *= static_cast<double>((nn(rp) + 1)) / nn(z(s)); 
            // }

            // update z(j)
            z(j) = sample_index(safe_exp(log_prob)); 
            // Rcout << "zj_new = " << z(j) << "\n";

            // update m and mbar tensors
            if (z(j) != zj_old) {
                m.slice(z(j)) += D;
                mbar.slice(z(j)) += Dbar;
                m.slice(zj_old) -= D;
                mbar.slice(zj_old) -= Dbar;
            }
        }
        
        // record label history
        z_hist.col(iter) = z + 1; // +1 is to put the labels on 1-based indexing
    }

    return z_hist;
    // return Rcpp::List::create( 
    //     Rcpp::Named("z") = z,
    //     Rcpp::Named("xi") = xi
    // );
}

// [[Rcpp::export]]
void comp_count_tensors(List A,
                std::vector<arma::uvec> xi, const int L, 
                arma::uvec z, 
                arma::cube &m, arma::cube &mbar) {

    for (int j=0; j < A.length(); j++) {
        List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
        arma::mat lambda = out["lambda"];
        arma::umat NN = out["NN"]; 
        m.slice(z(j)) += lambda;
        mbar.slice(z(j)) += NN - lambda;
    }

}

// [[Rcpp::export]]
arma::umat nsbm_z_update_cpp(List A,
                std::vector<arma::uvec> xi,
                arma::uvec z_init, 
                double alpha_eta = 1, double beta_eta = 5,
                double pi0 = 1, double w0 = 1,
                int niter = 50, int K = 5, int L = 10, bool verb = true,
                const double perturb = 1e-11,
                const bool rand_init = true, const bool seq_g_update = true) {

    int J = A.length();
    arma::uvec n(J);

    for (int j=0; j < J; j++) {
        n(j) = Rcpp::as<arma::sp_mat>(A[j]).n_rows;
    }

    // xi(j)[i] is in 0,1,...,L-1 where i ranges in 0,1,...,n(j)-1
    arma::uvec z(J); // z(k) is in 0,1,...,K-1
    arma::umat z_hist(J, J*niter, arma::fill::zeros);
    arma::vec  pi(K);
    arma::mat  w(L, K);
    // arma::cube beta_tensor(L, L, K, arma::fill::ones);

    arma::cube m(L, L, K, arma::fill::zeros); // m tensor
    arma::cube mbar(L, L, K, arma::fill::zeros); // mbar tensor

    z = z_init; 
    comp_count_tensors(A, xi, L, z, m, mbar); // updates m and mbar

    // updating z
    for (int iter = 0; iter < niter; iter++) {

        for (int j = 0; j < J; j++) {

            int zj_old = z(j);

            // calculate log_prob for updating z(j)
            List out = comp_blk_sums_and_sizes(A[j], xi[j], L);
            arma::mat D = out["lambda"]; 
            arma::mat Dbar = Rcpp::as<arma::umat>(out["NN"]) - D;

            arma::uvec xi_j_freq = get_freq(xi[j], L);
            arma::vec log_prob = comp_tensor_log_beta_ratio_sums(
                m, mbar, D, Dbar, zj_old, alpha_eta, beta_eta
            );
            // log_prob += log(w.t() + perturb) * xi_j_freq;
            // log_prob += log(pi + perturb);
        
            // update z(j)
            z(j) = sample_index(safe_exp(log_prob)); 

            // update m and mbar tensors
            if (z(j) != zj_old) {
                m.slice(z(j)) += D;
                m.slice(zj_old) -= D;
                mbar.slice(z(j)) += Dbar;
                mbar.slice(zj_old) -= Dbar;
            }

            z_hist.col(iter*J + j) = z + 1; 
        } // j
    }
    
    return z_hist;
    // return Rcpp::List::create( 
    //     Rcpp::Named("m") = m,
    //     Rcpp::Named("mbar") = mbar
    // );

  // sample_int_vec(K, J);
}