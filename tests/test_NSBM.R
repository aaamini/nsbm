source("R/data_gen.R")
source("R/inference.R")
source("R/competing_methods.R")
source("R/matching.R")
source("R/plotting.R")
source("R/nsbm_wrapper.R")
library(parallel)
library(ggplot2)
library(dplyr)
library(tibble)
library(Matrix)
library(RcppHungarian)
library(patchwork)

# Rcpp::sourceCpp("src/MCSBM.cpp", verbose = T)
Rcpp::sourceCpp("src/NSBM.cpp", verbose = T)
# Rcpp::sourceCpp("src/SBM.cpp", verbose = T)

# set.seed(1235)
n = 100
Ktru = 2
Ltru = 5
K = L = 20
J = 50

lambda = 20
# lambda = 30
zeta = .9 # try z = 0.3
gam = .7
nreps = 25
n_cores = 32
niter = 100

nathan_data = F

if (nathan_data) {  Ltru = 3 }

methods = list()

aa = 5
bb = 200
# aa = 1
# bb = 1
w0 = 0
pi0 = 1

methods[["NSBM"]] = function(A, K, L, niter=50) {
#    out = mix_mcsbm(A, K, L, 1, 1, 0.9, 1, niter, 3)
  mod = new(NSBM, A, K, L, aa, bb, w0, pi0) #, 0.9, 0.99)
  out = mod$run_gibbs(niter)
#   
   list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
}

# methods[["C++ (non-collapsed v2)"]] = function(A, K, L, niter=50) {
#   fit_nsbm(A, K, L, niter, collapsed = F, version = 2)
#   # model =  new(NestedSBM, A, K, L)
#   # model$run_gibbs_via_eta(niter) 
# }


# methods[["NSBM-rnd"]] = function(A, K, L, niter=50) {
#   #    out = mix_mcsbm(A, K, L, 1, 1, 0.9, 1, niter, 3)
#   mod = new(NSBM, A, K, L, aa, bb, w0, pi0) #, 0.9, 0.99)
#   mod$rnd_prob = 0.9
#   mod$decay = .99
#   out = mod$run_gibbs(niter)
#   #   
#   list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
# }

# methods[["NSBM-mix"]] =  function(A, K, L, niter=50)  {
#   out = mix_nsbm(A, K, L, 1, 1, 1, 1, decay = 0.95, sa_temp = 100, rnd_prob = 0.5, niter, 3)  
#   list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
# }

methods[["NSBM-SA"]] = # fit_mcsbm
 function(A, K, L, niter=50) {
#  mod = new(MCSBM, A, K, L, 1, 1, 0, 0.95)
  mod = new(NSBM, A, K, L, aa, bb, w0, pi0)
  mod$rnd_prob = 0
  mod$decay = 0.99
  mod$sa_temp = 100
  out = mod$run_gibbs(niter)

  list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
}

# methods[["MCSBM-tru"]] = function(A, K, L, xi_tru, niter=50) {
#   mod = new(MCSBM, A, K, L, 1, 1)
#   mod$set_xi_to_given(xi_tru)
#   mod$update_count_tensors()
#   
#   z_hist = matrix(0, nrow = J, ncol = niter+1)
#   z_hist[,1] = mod$z + 1
#   for (iter in 1:niter) {
#     mod$update_eta()
#     # round(mod$eta,2)
#     # simplify2array(eta)
#     
#     mod$update_w()
#     mod$update_pi()
#     for (j in 1:J) {
#       #for (s in 1:nrow(A[[j]])){
#       mod$update_z_element(j-1)
#       #}
#     }
#     z_hist[,iter+1] = mod$z + 1
#   }
#   list(z = z_hist, xi = lapply(1:(niter+1), function(dummry) xi_tru))
#   # out = mod$run_gibbs(niter)
#   # list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
# }


# out = gen_rand_nsbm(n=n, K=K, L=L, J=J,  lambda=lambda, gam=gam, zeta=zeta, sort_z = T)
# out = generate_nathans_data(n = n, J = J)
# A = out$A
# z_tru = out$z
# xi_tru = out$xi
# mout = methods[["MCSBM"]](A, K, L, 100)
# z_hist = mout$z
# xi_hist = mout$xi
# hsbm::get_slice_nmi(xi_hist[[51]], xi_tru)

mtd_names = names(methods)

# bench::mark(v1 = fit_sbm(A[1:3], L), v2 = fit_musbm(A[1:3], L), check = F)

res = do.call(rbind, mclapply(1:nreps, function(rep) {
# res = do.call(rbind, lapply(1:nreps, function(rep) {
  if (nathan_data) {
    # out = generate_nathans_data(n = n, J = J)
    # out = generate_nathans_data(n = n, J = J, lambda = lambda)
    out = generate_nathans_data(n = n, J = J, K = Ktru, lambda = lambda)  
  }
  else {
    # out = gen_rand_nsbm(n=n, K=Ktru, L=Ltru, J=J, lambda=lambda, gam=gam, zeta=zeta, sort_z = T)
    out = gen_rand_graphon(n = n, J = J, lambda = lambda, L = Ltru, std = 10)    
  }
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
    # eta = out$eta
  # mean_deg = mean(sapply(A, function(Adj) mean(rowSums(Adj))))
  
  do.call(rbind, lapply(seq_along(methods), function(j) { 
    # dt = as.numeric(system.time( mout <- methods[[j]](A, K, L, xi_tru, niter) )["elapsed"])
    dt = as.numeric(system.time( mout <- methods[[j]](A, K, L, niter) )["elapsed"])
    z_hist = mout$z
    xi_hist = mout$xi
    data.frame(method = mtd_names[j], 
               rep = rep,
               # mean_deg = mean_deg,
               iter = 1:(niter + 1), 
               z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru)),
               matching_score = sapply(xi_hist, function(xi) comp_matching_score(xi, xi_tru, z_tru)),
               xi_nmi = sapply(xi_hist, function(xi) hsbm::get_slice_nmi(xi, xi_tru)),
               elapsed_time = dt)
  }))
}, mc.cores = n_cores))
# }))    



res = as_tibble(res) %>% mutate(rep = as.character(rep))
# state_str =  sprintf("J = %d, n = %d, nr = %d, lam = %s, gam = %2.2f", 
#                      J, n, nreps, if(is.null(lambda)) "nathan" else lambda, gam)
state_str =  sprintf("J = %d, n = %d, nr = %d, lam = %s, K0 = %d, L0 = %d, K = %d, L = %d, (a,b) = (%2.1f, %2.1f)", 
                     J, n, nreps, if(is.null(lambda)) "NA" else lambda, Ktru, Ltru, K, L, aa, bb)

if (nathan_data) {
  state_str = sprintf("%s_nathan", state_str)
} else {
  state_str = sprintf("%s, zeta = %1.1f, gam = %1.1f", state_str, zeta, gam)
}

p1 = plot_paths_and_avg(res, xi_nmi) + ylab("xi-NMI") + labs(title = state_str) 
p2 = plot_paths_and_avg(res, z_nmi, alpha_range = c(0.2, 1)) + ylab("z-NMI") # + labs(title = state_str)
p3 = plot_paths_and_avg(res, matching_score) + ylab("Matching score") # + labs(title = state_str)

print(p1 + p2 + p3)

res %>% group_by(method) %>% summarise(avg_time = mean(elapsed_time))
# ggsave(sprintf("test_nsbm3_%s.png", state_str), width = 10, height=5)
