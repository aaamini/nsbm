// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "sampling.h"
#include "beta_calcs.h"

using namespace Rcpp;


void print_progress(int itr, int itr_max) {
  int width = ceil(log10(itr_max));
  if (itr % 10 == 0) Rcout << "("
                            << std::setw(width) << itr 
                            << " / " 
                            << std::setw(width) << itr_max << ")\r";
}


// void print_arma_obj(auto & x) {
    // Rcpp::print(Rcpp:wrap(x));
// }



//' @export
// [[Rcpp::export]]
double nhamming(const arma::sp_mat & A, int i, int j) {
  
    arma::sp_mat Aij = arma::join_rows(A.col(i), A.col(j));
    
    // arma::sp_mat::const_row_iterator it = A.begin_row();
    double ham = 0;
    for (int r = 0; r < A.n_rows; r++) {
        ham += Aij(r,0) != Aij(r,1);
    }
    //for (; it != A.end_row(); it++) {
    //    Rcout << (*it) << "\n";
    // }
    
    return ham / A.n_cols;
}

// [[Rcpp::export]]
arma::vec get_multi_nhamming(std::vector<arma::sp_mat> A, arma::umat index_list, 
                          const int type = 1) {
    int m = index_list.n_rows;
    int nlayers = A.size();
    arma::vec total_nham(m, arma::fill::zeros);
    for (int r = 0; r < m; r++) {
        arma::vec temp(nlayers, arma::fill::zeros);
        for (int t = 0; t < nlayers; t++) {
        //  arma::sp_umat At = A[t];
            temp[t] = nhamming(A[t], index_list(r,0), index_list(r,1));
            
         }
         switch (type) {
             case 0:
                total_nham[r] = temp.min();
                break;
             case 2:
                total_nham[r] = temp.max();
                break;
             default:
                total_nham[r] = arma::mean(temp);                
         }
         
    }
    return total_nham;
}

// List gem_posterior_counts(arma::uvec z, int K) {
//     arma::uvec count1(K, arma::fill::zeros);
//     arma::uvec count2(K, arma::fill::zeros);
//     int n = z.n_elem;

    
//     for (int k = 0; k < K; k++){
//       for (int i = 0; i < n; i++){
//         if (z(i) == k) {
//           count1(k)++;
//         } else if (z(i) > k) {
//           count2(k)++;
//         }
//       } //i
//     }      

//     return List::create(Named("count1") = count1, Named("count2") = count2);
// }

// [[Rcpp::export]]
arma::uvec get_up_freq(arma::uvec freq) {
    
    const int K = freq.n_elem;
    arma::uvec up_freq(K, arma::fill::zeros);

    for (int k = K-2; k >= 0; k--) {
        up_freq(k) = up_freq(k+1) + freq(k+1);
    }
    return up_freq;
}


// [[Rcpp::export]]
int find_tunc(arma::vec beta, double threshold) {
  int n = beta.n_elem;
  int idx = n;
  double cumsum = 0.0;
  for (idx = 0; idx < n - 1; idx++) {
    cumsum += beta(idx);
    if (cumsum > 1 - threshold) break;
  }
  return idx;
}


// [[Rcpp::export]]
arma::vec fast_agg(arma::vec x, arma::uvec z, int K) {
    int n = z.n_elem;

    arma::vec S(K, arma::fill::zeros);

    for (int i = 0; i < n; i++) {
        S(z(i)) += x(i);
    }
    return S;
}

arma::uvec fast_agg_u(arma::uvec x, arma::uvec z, int K) {
    int n = z.n_elem;

    arma::uvec S(K, arma::fill::zeros);

    for (int i = 0; i < n; i++) {
        S(z(i)) += x(i);
    }
    return S;
}


// counts the occurence frequency of integers 0, 1, 2, ..., K-1 in z
// the numbers in z should be in [0, K-1], otherwise an error occurs
// [[Rcpp::export]]
arma::uvec get_freq(arma::uvec z, int K) { 
  int n = z.n_elem;
  // int K = max(z)+1;
  arma::uvec freq(K,  arma::fill::zeros);
  // Rcout << K << freq; 
  for (int i = 0; i < n; i++) {
    freq(z(i))++;
  }
  return freq;
}

// [[Rcpp::export]]
arma::umat get_freq_minus_self(arma::uvec z, int K) {
    int n = z.n_elem;
    arma::uvec zcounts = get_freq(z, K);
    arma::umat freq = arma::repmat(zcounts.t(), n, 1);

    for (int i = 0; i < n; i++) {
        freq(i, z(i))--;
    }
    return freq;
}


// [[Rcpp::export]]
arma::mat comp_blk_sums(arma::sp_mat At, arma::uvec z, int Kcap) {
    // Compute block sums of a sparse matrix At w.r.t. labels in "z". The labels
    // in z are in the interval [0,1,...,Kcap]

    arma::sp_mat::const_iterator it     = At.begin();
    arma::sp_mat::const_iterator it_end = At.end();

    arma::mat lambda(Kcap, Kcap, arma::fill::zeros); 
    for(; it != it_end; ++it) {
        lambda(z(it.row()), z(it.col())) += (*it);
    }

    return lambda;
}

// [[Rcpp::export]]
arma::mat sp_compress_col(arma::sp_mat At, arma::uvec z, int Kcap) {
    int n = At.n_rows;
    arma::mat B(n, Kcap, arma::fill::zeros);

    arma::sp_mat::const_iterator it     = At.begin();
    arma::sp_mat::const_iterator it_end = At.end();

    for(; it != it_end; ++it) {
        B(it.row(), z(it.col())) += (*it);
    }

    return B;
}

// [[Rcpp::export]]
arma::vec sp_single_col_compress(arma::sp_mat A, int col_idx, arma::uvec z, int Kcap) {

    arma::sp_mat Acol(A.col(col_idx));
    arma::vec b(Kcap, arma::fill::zeros);

    // Rcpp::print(wrap(Acol));
    for (arma::sp_mat::iterator it = Acol.begin(); it != Acol.end(); ++it) {
        b(z(it.row())) += (*it);
        //Rcout << it.row() << " ";
    }
    return b;
}

// // This is for DCSBM setup
// arma::mat comp_blk_sums_diff(arma::sp_mat& A, int s, int zs_new, arma::uvec& z, int Kcap) {
//     // This requires fixing since it is assumming the Poi-DCSBM setup 

//     int zs = z(s);
//     double Ass =  A(s,s);

//     arma::vec U =  sp_single_col_compress(A, s, z, Kcap);
//     U(z(s)) -= Ass;

//     arma::vec delta(Kcap, arma::fill::zeros);
//     delta(zs_new)++;
//     delta(zs)--;

//     arma::mat D = delta*U.t() + U*delta.t();

//     // We can put an if here: if (zs != zs_new) -- not sure if it improves performance
//     D(zs_new, zs_new) += Ass;
//     D(zs, zs) -= Ass;

//     return D;
// }

// [[Rcpp::export]]
arma::mat comp_blk_sums_diff_v1(const arma::vec& U, const int zs_new, const int zs_old) {
    // Assumes A(s,s) is zero
    // arma::vec U =  sp_single_col_compress(A, s, z, Kcap);

    arma::vec delta(U.n_elem, arma::fill::zeros);
    delta(zs_new)++;
    delta(zs_old)--;

    arma::mat D = delta*U.t() + U*delta.t();
    D.diag() /= 2;
    
    return D;
}


// [[Rcpp::export]]
arma::mat comp_blk_sums_diff_v2(const arma::vec& U, const int zs_new, const int zs_old) {
    // Assumes A(s,s) is zero
    // arma::vec U = sp_single_col_compress(A, s, z, K);

    arma::mat D(U.n_elem, U.n_elem, arma::fill::zeros);
    D.row(zs_new) = U.t();
    D.row(zs_old) = -U.t();
    D.col(zs_new) += U;
    D.col(zs_old) -= U;
    D.diag() /= 2;
    
    return D;
}


template<typename Func>
double sym_prod(Func f, const int rp, const int r, const int K) {
    // f is a symmetric function
    double prod_result = 1;
    for (int y = 0; y < K; y++) {
        prod_result *= f(r,y) * f(rp,y);
    }
    double temp = f(r,rp);
    if (temp > 0) {
        prod_result /= temp;
    }
    return prod_result;
    // return prod_result / f(r,rp);
}

template<typename Func>
double sym_sum(Func f, const int rp, const int r, const int K) {
    // f is a symmetric function
    double result = 0;
    for (int y = 0; y < K; y++) {
        result += f(r,y) + f(rp,y);
    }
    return result - f(r,rp);;
}


// double comp_beta_ratio_prod(
//     const arma::mat& m, 
//     const arma::mat& mbar, 
//     const arma::vec& U,
//     const arma::uvec& V, 
//     const int zs_new, const int zs_old,
//     const int alpha, const int beta) {

//     arma::mat D = comp_blk_sums_diff_v2(U, zs_new, zs_old);
//     arma::mat DN = comp_blk_sums_diff_v2(arma::conv_to<arma::vec>::from(V), zs_new, zs_old);
//     arma::mat m_new = m + D;
//     arma::mat mbar_new = mbar + DN - D;

//     // A lambda, requires C++11
//     auto f = [&](int x, int y) {
//         return 
//             R::beta(m_new(x, y) + alpha, mbar_new(x, y) + beta) / 
//             R::beta(m(x, y) + alpha, mbar(x, y) + beta);
//     };
//     // Rcpp::Rcout << f(1,2);

//     return sym_prod(f, zs_new, zs_old, m.n_cols);
// }


// [[Rcpp::export]]
arma::vec comp_beta_ratio_prods_v1(
    const arma::mat& m, 
    const arma::mat& mbar, 
    const arma::vec& U,
    const arma::uvec& V, 
    const int zs_old,
    const int alpha, const int beta) {

    int K = m.n_cols;
    arma::vec out(K);

    for (int zs_new = 0; zs_new < K; zs_new++) {

        if (zs_new == zs_old) {
            out[zs_new] = 1;
            continue;
        }
        arma::mat D = comp_blk_sums_diff_v1(U, zs_new, zs_old);
        arma::mat DN = comp_blk_sums_diff_v1(arma::conv_to<arma::vec>::from(V), zs_new, zs_old);
        arma::mat m_new = m + D;
        arma::mat mbar_new = mbar + DN - D;
        arma::mat Dbar = DN - D;

        // print(wrap(m));
        // print(wrap(m_new));
        // print(wrap(mbar));
        // print(wrap(mbar_new));

        // A lambda, requires C++11
        auto f = [&](int x, int y) {
            return
                comp_beta_ratio_v1(m(x,y) + alpha, mbar(x,y) + beta, D(x,y), Dbar(x,y));
                // R::beta(m_new(x, y) + alpha, mbar_new(x, y) + beta) / 
                //  ( R::beta(m(x, y) + alpha, mbar(x, y) + beta) + DBL_MIN );
        };
        // Rcpp::Rcout << f(1,2);

        out(zs_new) = sym_prod(f, zs_new, zs_old, K);

    } // zs_new
    
    return out;
}

// [[Rcpp::export]]
arma::vec comp_log_beta_ratio_sums(
    const arma::mat& m, 
    const arma::mat& mbar, 
    const arma::vec& U,
    const arma::uvec& V, 
    const int zs_old,
    const int alpha, const int beta) {

    int K = m.n_cols;
    arma::vec out(K);

    for (int zs_new = 0; zs_new < K; zs_new++) {

        if (zs_new == zs_old) {
            out[zs_new] = 0;
            continue;
        }
        arma::mat D = comp_blk_sums_diff_v1(U, zs_new, zs_old);
        arma::mat DN = comp_blk_sums_diff_v1(arma::conv_to<arma::vec>::from(V), zs_new, zs_old);
        // arma::mat m_new = m + D;
        // arma::mat mbar_new = mbar + DN - D;
        arma::mat Dbar = DN - D;

        // print(wrap(m));
        // print(wrap(m_new));
        // print(wrap(mbar));
        // print(wrap(mbar_new));

        // A lambda, requires C++11
        auto f = [&](int x, int y) {
            return
                comp_log_beta_ratio(m(x,y) + alpha, mbar(x,y) + beta, D(x,y), Dbar(x,y));
                // R::beta(m_new(x, y) + alpha, mbar_new(x, y) + beta) / 
                //  ( R::beta(m(x, y) + alpha, mbar(x, y) + beta) + DBL_MIN );
        };
        // Rcpp::Rcout << f(1,2);

        out(zs_new) = sym_sum(f, zs_new, zs_old, K);

    } // zs_new
    
    return out;
}


// This is just a test function -- to be removed


// arma::vec test_N_update(arma::sp_mat& A, arma::mat N, int s, arma::uvec z, int Kcap, double alpha, double beta) {
//     arma::vec prob(Kcap, arma::fill::zeros);
//     arma::mat Np(Kcap, Kcap, arma::fill::zeros);

//     arma::uvec nn = get_freq(z, Kcap);
    
//     for (int rp = 0; rp < Kcap; rp++) {
//         // Rcpp::Rcout << update_blk_sums(At, N, s, rp, z, Kcap);
//         // Np = update_blk_sums(At, N, s, rp, z, Kcap);
//         Np = N + comp_blk_sums_diff(A, s, rp, z, Kcap);
//         prob(rp) = ratio_fun(N.as_col(), Np.as_col(), alpha, beta);
//         //Rcpp::Rcout << (nn(rp) + 1) / nn(z(s));
//         prob(rp) *= static_cast<double>((nn(rp) + 1)) / nn(z(s)); 
//     }
//     return prob;
// }



// arma::vec prod_dist(arma::sp_mat At, arma::mat N, int s, arma::uvec z, int Kcap, double alpha, double beta) {
//     arma::vec out(Kcap, arma::fill::zeros);
//     arma::mat Np(Kcap, Kcap, arma::fill::zeros);

    
//     for (int rp = 0; rp < Kcap; rp++) {
//         Np = update_blk_sums(At, N, s, rp, z, Kcap);
//         out(rp) = ratio_fun(N.as_col(), Np.as_col(), alpha, beta);
//     }
//     return out;
// }


// TODO: This should be named comp_blk_sums_and_sizes
// [[Rcpp::export]]
List comp_blk_sums_and_sizes(arma::sp_mat At, arma::uvec z, int Kcap, bool div_diag = true) {
    // z is a Kcap x 1 vector
    // At is a sparse n x n matrix

    int n = At.n_rows;
    // arma::mat lambda(Kcap, Kcap, arma::fill::zeros);

    arma::uvec zc = get_freq(z, Kcap); // zcounts

    arma::mat lambda = comp_blk_sums(At, z, Kcap);
    arma::umat NN = zc * zc.t() - arma::diagmat(zc);

    if (div_diag) {
        lambda.diag() /= 2; 
        NN.diag() /= 2;     // assumes that the diagonal of At is zero, 
                            // otherwise have to remove diag. first
    }
    
    return List::create( 
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("NN") = NN
    );
}


std::pair<arma::mat, arma::mat> comp_blk_sums_and_comps(
    const arma::sp_mat& At, const arma::uvec &z, int Kcap, bool div_diag = true) {
    // z is a Kcap x 1 vector
    // At is a sparse n x n matrix

    int n = At.n_rows;
    // arma::mat lambda(Kcap, Kcap, arma::fill::zeros);

    arma::uvec zc = get_freq(z, Kcap); // zcounts

    arma::mat lambda = comp_blk_sums(At, z, Kcap);
    arma::umat NN = zc * zc.t() - arma::diagmat(zc);

    if (div_diag) {
        lambda.diag() /= 2; 
        NN.diag() /= 2;     // assumes that the diagonal of At is zero, 
                            // otherwise have to remove diag. first
    }
    
    return std::make_pair(lambda, NN-lambda);
}

// // [[Rcpp::export]]
// arma::mat test_comp_blk_sums_and_comps(std::pair<arma::mat, arma::mat> a_pair) {
//     auto [x, y] = a_pair;

//     List::create(
//         Named("x") = x,
//         Named["y"] = y
//     )
// }


// [[Rcpp::export]]
arma::mat beta_fun_symmat(arma::mat a, arma::mat b) {
    // Computes beta function over symmetric matrix arguments
    // only computes the lower triangular part + diagonal; set the rest to 1
    int K = a.n_rows;
    arma::mat out(K, K, arma::fill::ones);
    for (int k = 0; k < K; k++) {
        for (int el = 0; el <= k; el++) {
            out(k, el) = R::beta(a(k, el), b(k, el));
        }
    }
    return out;
}

// [[Rcpp::export]]
arma::mat comp_beta_matrix(const arma::sp_mat& A, arma::uvec& z, const int K, double alpha, double beta) {
    List out = comp_blk_sums_and_sizes(A, z, K);
    arma::mat lambda = out["lambda"];
    arma::umat NN = out["NN"]; 
            
    return beta_fun_symmat(lambda + alpha, NN - lambda + beta);
}




// An illustration of how to access only the zero elementv s of a sparse matrix
// [[Rcpp::export]]
void iter_over_sp_mat(arma::sp_mat At) {
    arma::sp_mat::const_iterator it     = At.begin();
    arma::sp_mat::const_iterator it_end = At.end();

    Rcout << std::setw(3) << "i"
        << std::setw(3) << "j"
        << std::setw(3) << "x"
        << std::endl;
    for(; it != it_end; ++it) {
        Rcout << std::setw(3) << it.row() + 1 
            << std::setw(3) << it.col() + 1
            << std::setw(3) << (*it) 
            << std::endl;
    }

}



