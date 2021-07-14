#set.seed(1234)
n = 50
K = 4
J = 6
z_tru = sample(1:K, J, replace = T)

L = 5
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = 10))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

source("nathans_nsbm/nSBM_functions.R")
out = gibbs.nSBM(A, ns = 10, monitor = T)
