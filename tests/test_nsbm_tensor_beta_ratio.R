# Clean test file -- should run without trouble

Rcpp::sourceCpp("src/utils.cpp", verbose = T)
#set.seed(1234)
n = 30
K = 4
L = 5
J = 10
alpha = 1
beta = 1

# data
z_tru = sample(1:K, J, replace = T)
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = 7))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

# tensor comp
comp_m_tensor = function(A, z, xi, K, L) {
    lapply(1:K, function(k) {
      idx = which(z == k)
      if (length(idx) == 0) return( matrix(0, L, L) )
      Reduce(`+`, lapply(idx, 
                         function(j) comp_blk_sums_and_sizes(A[[j]], xi[[j]]-1, L)$lambda))  
    }
  )
}

comp_mbar_tensor = function(A, z, xi, K, L) {
  lapply(1:K, function(k) {
    idx = which(z == k)
    if (length(idx) == 0) return( matrix(0, L, L) )
    Reduce(`+`, lapply(idx,  
                       function(j) {
                         out = comp_blk_sums_and_sizes(A[[j]], xi[[j]]-1, L)
                         out$NN - out$lambda
                       }))
    }
  )
}

upper_tri_sum = function(X) {
  sum(X - diag(diag(X)))/2 + sum(diag(X))
}

upper_tri_sum2 = function(X) {
  uidx = upper.tri(diag(nrow(X)), diag = T)
  sum(X[uidx])
}

## tests
m_tensor = comp_m_tensor(A, z_tru, xi_tru, K,  L)
mbar_tensor = comp_mbar_tensor(A, z_tru, xi_tru, K,  L)
q = simplify2array(m_tensor) # turn list into 3d array
qbar = simplify2array(mbar_tensor)

out = comp_blk_sums_and_sizes(A[[j]], xi_tru[[j]]-1, L)
D = out$lambda
Dbar = out$NN - D

log_prob_R = rep(0, K)
r0 = z_tru[j]
z = z_tru
for (r in 1:K) {
  z[j] = r
  q_new = simplify2array(comp_m_tensor(A, z, xi_tru, K,  L))
  qbar_new = simplify2array(comp_mbar_tensor(A, z, xi_tru, K,  L))
  
  temp = log(beta(q_new + alpha, qbar_new + beta)) - 
    log(beta(q + alpha, qbar + beta))
  # temp = log(beta(q_new+ alpha, qbar_new + beta))
  temp_k_sum = rowSums(temp, dims=2)
  log_prob_R[r] = upper_tri_sum2(temp_k_sum) 
}

Rcpp::sourceCpp("src/utils.cpp", verbose = T)

log_prob_cpp = as.vector(
  comp_tensor_log_beta_ratio_sums(q, qbar, D, Dbar, r0-1, alpha, beta)
)

#remove_max = function(x) x - max(x)
#rbind(remove_max(log_prob_R), remove_max(log_prob_cpp))


rbind(log_prob_R, log_prob_cpp)
cat(sprintf('Error = %e', sum(abs(log_prob_R - log_prob_cpp))))
