# set.seed(1234)
n = 100
K = 5
z = sample(1:K, n, replace = T)

eta = nett::gen_rand_conn(n, K, lambda = 10)
A = nett::fast_sbm(z, eta)


Rcpp::sourceCpp("src/utils.cpp", verbose = T)
 
out = comp_blk_sums_and_sizes(A, z-1, K)
out$lambda
z_new = z
s = sample(n, 1)
z_new[s] = sample(setdiff(1:K, z[s]), 1)
c(z[s], z_new[s])
out2 = comp_blk_sums_and_sizes(A, z_new-1, K)

U =  sp_single_col_compress(A, s-1, z-1, K);
(D1 = out2$lambda - out$lambda)
# (D2 = comp_blk_sums_diff_v1(A, s-1, z_new[s]-1, z-1, K))
# (D3 = comp_blk_sums_diff_v2(A, s-1, z_new[s]-1, z-1, K))
(D2 = comp_blk_sums_diff_v1(U, z_new[s]-1, z[s]-1))
(D3 = comp_blk_sums_diff_v2(U, z_new[s]-1, z[s]-1))


D1 - D2
D1 - D3


(DN1 = out2$NN - out$NN)
V = get_freq(z-1, K)
V[z[s]] = V[z[s]] - 1
(DN2 = comp_blk_sums_diff_v1(V, z_new[s]-1, z[s]-1))
(DN3 = comp_blk_sums_diff_v2(V, z_new[s]-1, z[s]-1))

DN1 - DN2
DN1 - DN3

alpha = 2
beta = 1
BetProd1 = prod(comp_beta_matrix(A, z_new-1, K, alpha, beta)/comp_beta_matrix(A, z-1, K, alpha, beta))

r = z[s]
rp = z_new[s]
idx = c(r,rp)
m = out$lambda
mbar = out$NN - m

m_new = m + D3
mbar_new = mbar + DN3 - D3

temp = beta(m_new + alpha, mbar_new + beta) / beta(m + alpha, mbar + beta)
BetProd2 = prod(temp[idx, ]) / prod(temp[r,rp])

BetProd1
BetProd2
# comp_beta_ratio_prod(m, mbar, U, V, rp-1, r-1, alpha, beta)
comp_beta_ratio_prods(m, mbar, U, V,  r-1, alpha, beta)


# out2$NN - out2$lambda - mbar_new

#bench::press(
  #iter = 1:20,
{
  n = 500
  K = 30
  z = sample(1:K, n, replace = T)
  
  eta = nett::gen_rand_conn(n, K, lambda = 10)
  A = nett::fast_sbm(z, eta)
  z_new = z
  s = sample(n, 1)
  z_new[s] = sample(setdiff(1:K, z[s]), 1)
  U =  sp_single_col_compress(A, s-1, z-1, K);
  
  bench::mark(
    min_iterations = 100,
    v0 = 
      comp_blk_sums_and_sizes(A, z_new-1, K)$lambda - 
      comp_blk_sums_and_sizes(A, z-1, K)$lambda,
    v1 = comp_blk_sums_diff_v1(U, z_new[s]-1, z[s]-1),
    v2 = comp_blk_sums_diff_v2(U, z_new[s]-1, z[s]-1),
    check = T
  )
}
#)

microbenchmark::microbenchmark(
  v0 = 
    comp_blk_sums_and_sizes(A, z_new-1, K)$lambda - 
    comp_blk_sums_and_sizes(A, z-1, K)$lambda,
  v1 = comp_blk_sums_diff_v1(U, z_new[s]-1, z[s]-1),
  v2 = comp_blk_sums_diff_v2(U, z_new[s]-1, z[s]-1),
  setup = {
    n = 500
    K = 12
    z = sample(1:K, n, replace = T)
    
    eta = nett::gen_rand_conn(n, K, lambda = 10)
    A = nett::fast_sbm(z, eta)
    z_new = z
    s = sample(n, 1)
    z_new[s] = sample(setdiff(1:K, z[s]), 1)
    U =  sp_single_col_compress(A, s-1, z-1, K)
  }
)
