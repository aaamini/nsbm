#set.seed(1234)
n = 50
K = 4
J = 10
z_tru = sample(1:K, J, replace = T)

L = 5
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = 10))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

Rcpp::sourceCpp("hsbm_package/src/utils.cpp", verbose = T)

# `comp_blk_sums` is a function in utils.cpp
comp_blk_sums(A[[2]], xi_tru[[2]]-1, L) + comp_blk_sums(A[[3]], xi_tru[[3]]-1, L)

comp_m_tensor = function(A, z, xi, L) {
  K = max(z)
  lapply(1:K, function(k) 
    Reduce(`+`, lapply(which(z == k), 
                       function(j) comp_blk_sums(A[[j]], xi[[j]]-1, L)))
  )
}

# Compute the true m tensor
m_tensor = comp_m_tensor(A, z_tru, xi_tru, L)
m_tensor

# Effect of changing a single xi
xi = xi_tru
s = 5
j = 8
xi[[j]][s] = (xi[[j]][s]  + 2) %% L + 1  # change a single label

c(xi[[j]][s], xi_tru[[j]][s] )

# # check that only one element is different
# lapply(1:J, function(j) sum(xi_tru[[j]] != xi[[j]]))

m_tensor_new = comp_m_tensor(A, z_tru, xi, L)
lapply(1:K, function(k) m_tensor_new[[k]] - m_tensor[[k]])

# The effect of changing a single z coordinate
z = z_tru
z[j] = (z_tru[j] + 1) %% K + 1
c(z[j], z_tru[j])
m_tensor_new = comp_m_tensor(A, z, xi_tru, L)
lapply(1:K, function(k) m_tensor_new[[k]] - m_tensor[[k]])
