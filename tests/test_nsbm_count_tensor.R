#set.seed(1234)
n = 30
K = 4
J = 10

z_tru = sample(1:K, J, replace = T)


L = 5
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = 7))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

Rcpp::sourceCpp("src/utils.cpp", verbose = T)

# `comp_blk_sums` is a function in utils.cpp
comp_blk_sums_and_sizes(A[[2]], xi_tru[[2]]-1, L)$lambda + comp_blk_sums_and_sizes(A[[3]], xi_tru[[3]]-1, L)$lambda

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

# Compute the true m tensor
m_tensor = comp_m_tensor(A, z_tru, xi_tru, K,  L)
mbar_tensor = comp_mbar_tensor(A, z_tru, xi_tru, K,  L)
m_tensor
mbar_tensor

# Effect of changing a single xi
xi = xi_tru
s = 5
j = 8
xi[[j]][s] = (xi[[j]][s]  + 2) %% L + 1  # change a single label

c(xi[[j]][s], xi_tru[[j]][s] )

# # check that only one element is different
# lapply(1:J, function(j) sum(xi_tru[[j]] != xi[[j]]))

m_tensor_new = comp_m_tensor(A, z_tru, xi, K, L)
lapply(1:K, function(k) m_tensor_new[[k]] - m_tensor[[k]])

# The effect of changing a single z coordinate
z = z_tru
# z[j] = (z_tru[j] + 1) %% K + 1
z[j] = 1
c(z[j], z_tru[j])


m_tensor_new = comp_m_tensor(A, z, xi_tru, K,  L)
lapply(1:K, function(k) m_tensor_new[[k]] / m_tensor[[k]])[[z_tru[j]]]

lapply(1:K, function(k) m_tensor_new[[k]] - m_tensor[[k]])




comp_blk_sums_and_sizes(A[[j]], xi_tru[[j]]-1, L)$lambda


##
alpha = 1
beta = 1


q = simplify2array(m_tensor)
qbar = simplify2array(mbar_tensor)

out = comp_blk_sums_and_sizes(A[[j]], xi_tru[[j]]-1, L)
D = out$lambda
Dbar = out$NN - D

upper_tri_sum = function(X) {
  sum(X - diag(diag(X)))/2 + sum(diag(X))
}

upper_tri_sum2 = function(X) {
  idx = upper.tri(diag(nrow(X)), diag = T)
  sum(X[uidx])
}


log_prob = log_prob2 = rep(0, K)
r0 = z_tru[j]
z = z_tru
for (r in 1:K) {
  z[j] = r
  q_new = simplify2array(comp_m_tensor(A, z, xi_tru, K,  L))
  qbar_new = simplify2array(comp_mbar_tensor(A, z, xi_tru, K,  L))
  
  temp = log(beta(q_new + alpha, qbar_new + beta)) - 
    log(beta(q + alpha, qbar + beta))
  temp_k_sum = rowSums(temp, dims=2)
  # temp = log(beta(m_tensor_new[[k]] + alpha, mbar_tensor_new[[k]] + beta))
  log_prob[r] = upper_tri_sum2(temp_k_sum) 
}

Rcpp::sourceCpp("src/utils.cpp", verbose = T)

log_prob2 = comp_tensor_log_beta_ratio_sums(q, qbar, D, Dbar, r0-1, alpha, beta)
# log_prob - max(log_prob)
# as.vector(log_prob2 - max(log_prob2))

log_prob
as.vector(log_prob2)
log_prob - as.vector(log_prob2)


# for (rp in 1:K) {
#   z[j] = rp
#   m_tensor_new = comp_m_tensor(A, z, xi_tru, K,  L)
#   mbar_tensor_new = comp_mbar_tensor(A, z, xi_tru, K,  L)
#   log_prob[rp] = sum(log(beta(m_tensor_new[[rp]] + alpha, mbar_tensor_new[[rp]] + beta)[uidx]))
#   
#   if (rp != r) {
#     log_prob2[rp] = sum(log(beta(m_tensor[[rp]] + D + alpha, mbar_tensor[[rp]] + DN + beta)[uidx]))
#   } else {
#     log_prob2[rp] = sum(log(beta(m_tensor[[rp]]+ alpha, mbar_tensor[[rp]] + beta)[uidx]))
#   }
#   
# }
