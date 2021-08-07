library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = F)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

#set.seed(1234)
n = 30
K = 4
L = 5
J = 20
alpha = 1
beta = 1
niter = 10
lambda = 7

# data
z_tru = sample(1:K, J, replace = T)
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = lambda))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

z_init = sample(1:K, J, replace = T) #
z_init_cpp = z_init - 1
xi_tru_cpp = lapply(xi_tru, function(x) x-1)

model = new(NestedSBM, A, K, L)
model

model$z
model$m
model$comp_count_tensors()
# model$xi[[10]]
model$m
model$A[[1]]
model$z
model$xi
model$w
model$pi

model$get_xi_freq_over_z(2)
model$update_w()
model$w
model$update_pi()
model$pi

model$xi = xi_tru_cpp
model$comp_count_tensors()
(m_old = model$m)
model$update_z_element(3)
model$m - m_old



## Set z to true z and update xi only
model$z = z_tru-1
model$z 
model$set_xi_to_random_labels()
xi_init = lapply(model$xi, function(x) x + 1)
xi_init
model$comp_count_tensors()
model$m

# Test whether "m" changes by xi-update
j = 5
(old_xi = model$xi[[j]])
old_m = model$m
model$update_xi_element(j-1, 3)
model$m - old_m
# model$xi[[1]] - old_xi

# Compate incremental update and full update
j = 8
old_m = model$m
model$update_xi_element(j-1, 3) # updates "m" by incremental update
m2 = model$m # updated m
m2 - old_m
model$comp_count_tensors() # full computation of the tensor
model$m - m2 # this should always be zero

model$blk_compressions
model$update_xi_element(j-1, 3)
out = model$run_gibbs(10)
