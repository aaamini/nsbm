source("R/data_gen.R")
source("R/inference.R")
source("R/competing_methods.R")
source("R/nsbm_wrapper.R")

sourceCpp("src/dpsbm.cpp")
library(parallel)
library(ggplot2)
library(dplyr)

set.seed(1234)
n = 400
K = 2
L = 5
J = 20
lambda = 30
zeta = .3
gam = .1
nreps = 10
n_cores = 32
niter = 100

out = gen_rand_nsbm(n=n, K=K, L=L, J=J,  lambda=lambda, gam = gam, zeta=zeta)
A = out$A
z_tru = out$z
xi_tru = out$xi

out = fit_nsbm(A, 15, 15, niter, collapsed = F)
z = get_map_labels(out$z)$labels
xi = get_map_labels(out$xi)$labels
nett::compute_mutual_info(z, z_tru)
hsbm::get_slice_nmi(xi, xi_tru)
hsbm::seq_nmi_plot(out$xi)
# 
out = spec_net_clust(A, K = K, L = L)
nett::compute_mutual_info(out$z, z_tru)
hsbm::get_slice_nmi(out$xi, xi_tru)
# 
xi = lapply(1:J, function(j) get_map_labels( fit_dpsbm(A[[j]], Zcap = L) )$label)
hsbm::get_slice_nmi(xi, xi_tru)

model =  new(NestedSBM, A, J, L)
model$z = 0:(J-1)
model$z
# model$w0 = 1
# model$pi0 = 1
colSums(model$w)
out = model$run_gibbs_via_eta(100)
out$xi  = lapply(out$xi, function(xi) lapply(xi, function(x) x + 1))
xi = get_map_labels(out$xi)$labels
z = get_map_labels(out$z)$labels
# xi
hsbm::get_slice_nmi(xi, xi_tru)

#model$set_beta_params(1, 1)
#model$w0 = 0.001
#model$pi0 = 0.001
# model

# nett::compute_mutual_info(ztru, sp_out$z)
# # get_map_labels(out$z)
# # get_map_labels(out$xi)$labels
# # 
# mean(Matrix::rowSums(A[[1]]))
# 
# K = L = 10
# niter = 50
# out = fit_nsbm(A, niter = niter, collapsed = F)
# 
# # library(nett)
# # out$z[,niter+1]
# # model =  new(NestedSBM, A, K, L)
# # fitted_model <- model$run_gibbs(niter) 
# # nett::compute_mutual_info(ztru, fitted_model$z[,niter+1])
# zh = out$z
# plot(sapply(1:ncol(zh), function(itr) nett::compute_mutual_info(ztru, zh[,itr])))
# 
# compute_mutual_info(ztru, get_map_labels(out$z)$labels)
# 
# j = 1
# out$xi[[1]]
# 
# sourceCpp("src/dpsbm.cpp")
# aa = fit_dpsbm(A[[1]], Zcap = 2*L)
# plot(sapply(1:ncol(aa), function(itr) nett::compute_mutual_info(xi[[1]], aa[,itr])))
# compute_mutual_info(xi[[1]], get_map_labels(aa)$labels)
# 
# compute_mutual_info(xi[[1]], spec_clust(A[[1]], 7))
# 
# # samp <- splice_sampler(A, K, L, ns = niter, monitor = TRUE)
# # samp$z


