n = 50
K = 5
z = sample(1:K, n, replace = T)
B = nett::gen_rand_conn(n, K, lambda = 10)
A = nett::fast_sbm(z, B)
diag(A) = sample(0:1, n, replace =T )  


# Rcpp::sourceCpp("hsbm_package/src/.cpp", verbose = T)
Rcpp::sourceCpp("hsbm_package/src/utils.cpp", verbose = T)
# x = runif(n)
# fast_agg(x, z-1, K)

s = 5
zp  = z 
zp[s] = (z[s] + 3) %% K + 1
A[s,s]
c(z[s],zp[s])

u = get_freq(z-1, K) 
M = u %*% t(u)
up = get_freq(zp-1, K) 
Mp = up %*% t(up)
M - Mp


library(Matrix)
  
(N = comp_blk_sums(A, z-1, K))
(Np = comp_blk_sums(A, zp-1, K))
(D1 = Np - N)
(D2 = nett::compute_block_sums(A, zp) - nett::compute_block_sums(A, z)) 

r = z[s]
rp = zp[s]
(U = fast_agg(A[s,], z-1, K))
U[z[s]] = U[z[s]] - A[s,s] 

delta = matrix(nett::label_vec2mat(zp[s],K) - nett::label_vec2mat(z[s],K), ncol=1) 
(D3 = U%*%t(delta) + delta %*% t(U))

(D4 = comp_blk_sums_diff(A, s-1, rp-1, z-1, K)) 

if (r != rp) {
  D3[rp, rp] =  D3[rp, rp] + A[s,s]
  D3[r, r] =  D3[r, r] - A[s,s]
}
  
D1-D2 
D2-D3
D3-D4

# Np
# N + D4

# zp[s] = 5
# (Np1 = comp_blk_sums(A, zp-1, K))
# N - Np1
# zp[s] = 3
# (Np2 = comp_blk_sums(A, zp-1, K))
# N - Np2


alpha1 = 0.5
beta1 = 0.25
fun = function(x) gamma(alpha1+x) / (beta1 +0.5)^(alpha1 + x)

u = fun(Np[1:5]/2)
v = f_fun(Np[1:5]/2, alpha1, beta1)
abs(u - v)/v



# microbenchmark::microbenchmark(a =  prod(fun(Np/2)/fun(N/2)), 
#                                b = prod(fun(Np[idx,-idx]/2)/fun(N[idx,-idx]/2))^2 * 
#                                  prod(fun(Np[idx,idx]/2)/fun(N[idx,idx]/2)))

a = prod(fun(Np/2)/fun(N/2))
a
ratio_fun(N, Np, alpha1, beta1)


Rcpp::sourceCpp("hsbm_package/src/utils.cpp", verbose = T)

# pd_cpp = t(prod_dist(A, N, s-1, z-1, K, alpha1, beta1))
# 
# pd_R = sapply(1:K, function(rp) {
#   zp  = z 
#   zp[s] = rp
#   Np = comp_blk_sums(A, zp-1, K)
#   prod(fun(Np/2)/fun(N/2))
# })
# 
# pd_cpp - pd_R

# pd2_cpp = t(test_N_update(A, N, s-1, z-1, K, alpha1, beta1))
# 
# nn = table(z)
# pd2_R = sapply(1:K, function(rp) {
#   zp  = z 
#   zp[s] = rp
#   Np = comp_blk_sums(A, zp-1, K)
#   prod(fun(Np/2)/fun(N/2)) * (nn[rp] + 1) / nn[z[s]]
# })
# pd2_R - pd2_cpp

Rcpp::sourceCpp("hsbm_package/src/collapsed_sbm.cpp", verbose = T)
As = A
out = fit_sbm(As, 5, niter = 10, alpha= 0.1, beta = 0.1)
t = 5
out$z_hist
table(out$z_hist[,t]+1)
out$nn_hist[,t]
out$N_hist[,t]
as.vector(nett::compute_block_sums(As, out$z_hist[,t]+1))


# (idx = c(r,rp))
# c = prod(fun(Np[idx,-idx]/2)/fun(N[idx,-idx]/2))^2 * 
#   prod(fun(Np[idx,idx]/2)/fun(N[idx,idx]/2))
# 
# a
# c
# 
# # b = (prod(fun(Np[r,-r]/2)/fun(N[r,-r]/2))*
# #        prod(fun(Np[-rp,rp]/2)/fun(N[-rp,rp]/2)))^2*
# #   fun(Np[r,r]/2)/fun(N[r,r]/2)*
# #   fun(Np[rp,rp]/2)/fun(N[rp,rp]/2)
