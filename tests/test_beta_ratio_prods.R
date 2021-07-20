n = 100
K = 5
z = sample(1:K, n, replace = T)

eta = nett::gen_rand_conn(n, K, lambda = 10)
A = nett::fast_sbm(z, eta)

Rcpp::sourceCpp("src/utils.cpp", verbose = T)
out = comp_blk_sums_and_sizes(A, z-1, K)
m = out$lambda
mbar = out$NN - m

s = sample(n, 1)
U =  sp_single_col_compress(A, s-1, z-1, K);
V = get_freq(z-1, K)
V[z[s]] = V[z[s]] - 1


alpha = 2
beta = 1
BetProd1 = BetProd2 = rep(0,K) 
for (zs_new in 1:K) {
  z_new = z
  z_new[s] = zs_new # sample(setdiff(1:K, z[s]), 1)
  BetProd1[zs_new] = prod(comp_beta_matrix(A, z_new-1, K, alpha, beta)/comp_beta_matrix(A, z-1, K, alpha, beta))
  # BetProd2[zs_new] = comp_beta_ratio_prod(m, mbar, U, V, zs_new-1, z[s]-1, alpha, beta)
}
(logBetSums1 = log(BetProd1))
(logBetSums3 = log(as.vector(comp_beta_ratio_prods_v1(m, mbar, U, V,  z[s]-1, alpha, beta))))
# as.vector(comp_beta_ratio_prods_v2(m, mbar, U, V,  z[s]-1, alpha, beta)) 
# abs(logBetSums3 - logBetSums1)/logBetSums1
logBetSums3 - logBetSums1
t(comp_log_beta_ratio_sums(m, mbar, U, V,  z[s]-1, alpha, beta))

# bench::mark(
#   min_iterations = 500,
#   v1 = comp_beta_ratio_prods(m, mbar, U, V,  z[s]-1, alpha, beta),
#   v2 = comp_beta_ratio_prods_v2(m, mbar, U, V,  z[s]-1, alpha, beta),
#   check = T
# )

##
Rcpp::sourceCpp("src/multsbm.cpp", verbose = T)
z_list = convert_cpp_label_matrix_to_list(
  multsbm_collapsed_gibbs_sampler_v2(A, K, alpha = alpha, beta = beta, niter = niter)
)
library(readr)
z = read_csv("z.csv", col_names = F)$X1
s = 79
z[s+1]
U = sp_single_col_compress(A, s, z, K)
V = get_freq(z, K)
V[z[s+1]+1] = V[z[s+1]+1] - 1 
U
V
read_csv("V.csv", col_names = F)$X1
read_csv("U.csv", col_names = F)$X1

(m = as.matrix(read_csv("m.csv", col_names = F)))
(mbar = as.matrix(read_csv("mbar.csv", col_names = F)))

out = comp_blk_sums_and_sizes(A, z, K)
out$lambda
out$NN - out$lambda 

Rcpp::sourceCpp("src/utils.cpp", verbose = T)
comp_beta_ratio_prods(m, mbar, U, V, z[s+1], alpha, beta)

res = rep(0,K)
for (rp in 0:2){
  znew = z
  znew[s+1] = rp
  out = comp_blk_sums_and_sizes(A, znew, K)
  m_new = out$lambda
  mbar_new = out$NN - out$lambda
  temp = beta(m_new + alpha, mbar_new + beta) / beta(m + alpha, mbar + beta)
  res[rp+1] = prod(temp[upper.tri(temp, diag = T)])
}
res
