source("R/data_gen.R")
source("R/inference.R")
source("R/competing_methods.R")
library(parallel)
library(ggplot2)
library(dplyr)
library(Matrix)
library(RcppHungarian)
library(patchwork)

Rcpp::sourceCpp("src/MuSBM.cpp", verbose = T)
Rcpp::sourceCpp("src/SBM.cpp", verbose = T)


# Caompute the optimal matching permutation from zh to z, both label vectors
matching_perm = function(zh, z) {
  n = length(z)
  out = HungarianSolver(1-compute_confusion_matrix(zh, z)/n)
  perm = out$pairs[,2]
  perm
}

comp_all_matching_perms = function(xih, xi_tru) {
  lapply(seq_along(xih), function(j) matching_perm(xih[[j]], xi_tru[[j]]))
}
comp_matching_score = function(xih, xi_tru){
  length(unique(comp_all_matching_perms(xih, xi_tru))) # / length(xih)
}



# set.seed(1235)
n = 200
K = 1 # K has to be 1 here
L = 5
J = 30
lambda = 30
zeta = .3 # This does not matter

gam = .7
nreps = 10
n_cores = 32
niter = 100



methods = list()
methods[["SBM"]] = function(A, L, niter=100) {
  temp_label_list = lapply(seq_along(A), function(j) {
    mod = new(SBM, A[[j]], L, 1, 1)
    mod$run_gibbs(niter) 
  })
  lapply(1:(niter+1), function(iter) 
    lapply(seq_along(temp_label_list), function(j) temp_label_list[[j]][,iter])
  )
}

methods[["MuSBM"]] = function(A, L, niter=100) {
  mod = new(MultSBM, A, L, 1, 1)
#   mod$set_xi_to_random_labels()
  temp_label_list = mod$run_gibbs(niter)
  lapply(temp_label_list, add_one_to_xi)
}

mtd_names = names(methods)

# bench::mark(v1 = fit_sbm(A[1:3], L), v2 = fit_musbm(A[1:3], L), check = F)

res = do.call(rbind, mclapply(1:nreps, function(rep) {
  # res = do.call(rbind, lapply(1:nreps, function(rep) {
  out = gen_rand_nsbm(n=n, K=K, L=L, J=J,  lambda=lambda, gam=gam, zeta=zeta, sort_z = T)
  A = out$A
  # z_tru = out$z
  xi_tru = out$xi
  eta = out$eta[[1]]
  
  do.call(rbind, lapply(seq_along(methods), function(j) { 
    dt = as.numeric(system.time( xi_hist <- methods[[j]](A, L, niter) )["elapsed"])
    
    data.frame(method = mtd_names[j], 
               rep = rep,
               iter = 1:(niter + 1), 
               matching_score = sapply(xi_hist, function(xi) comp_matching_score(xi, xi_tru)),
               xi_nmi =  sapply(xi_hist, function(xi) hsbm::get_slice_nmi(xi, xi_tru)),
               elapsed_time = dt)
  }))
}, mc.cores = n_cores))
# }))    

state_str =  sprintf("J = %d, n = %d, nr = %d, lam = %d, gam = %2.2f", J, n, nreps, lambda, gam)
p1 = res %>% 
  group_by(iter, method) %>% summarise(xi_nmi = mean(xi_nmi)) %>% 
  ggplot(aes(x = iter, y = xi_nmi, color = method)) + 
  # geom_line(aes(size = method), alpha = 0.5) +
  geom_line(size = 1.2) + 
  theme_minimal() +
  # scale_colour_manual(values = c(1,1.5)) +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("xi-NMI") + xlab("Iteration") + labs(title = state_str)

p2 = res %>% 
  group_by(iter, method) %>% summarise(matching_score = mean(matching_score)) %>% 
  ggplot(aes(x = iter, y = matching_score, color = method)) + 
  # geom_line(aes(size = method), alpha = 0.5) +
  geom_line(size = 1.2) + 
  theme_minimal() +
  # scale_colour_manual(values = c(1,1.5)) +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("Matching score") + xlab("Iteration") + labs(title = state_str)



p1 + p2
# ggsave(sprintf("test_musbm_%s.png", state_str), width = 10, height=5)

# out = gen_rand_nsbm(n=n, K=K, L=L, J=J,  lambda=lambda, gam = gam, zeta=zeta, sort_z = T)
# A = out$A
# z_tru = out$z
# xi_tru = out$xi
# eta = out$eta[[1]]

# # check whether eta is permutation invariant
# rand_perm = sample.int(L)
# rand_P = as(rand_perm, "pMatrix")
# sum((rand_P %*% eta %*% t(rand_P) - eta)^2)


# # xih_list1 = lapply(seq_along(A), function(j) fit_sbm(A[[j]], L, niter))
# xih_list1 = fit_sbm(A, L)
# # xih_list1 = lapply(1:(niter+1), function(iter) 
# #    lapply(seq_along(xih_list1), function(j) xih_list1[[j]][,iter])
# #  )
# 
# 
# xih_list2 = fit_musbm(A, L, niter)
# # xih_list2 = lapply(xih_list2, add_one_to_xi)

# compute_nmi_evolution = function(xih_list) {
#   sapply(seq_along(xih_list), function(iter) 
#     hsbm::get_slice_nmi(xih_list[[iter]], xi_tru))
# }
# 
# iter_seq = 1:(niter+1)
# res = data.frame(iter = iter_seq, nmi = compute_nmi_evolution(xih_list1), method = "SBM")
# res = rbind(res, 
#             data.frame(iter = iter_seq, nmi = compute_nmi_evolution(xih_list2), method = "MuSBM")
# )

# res %>% 
#   ggplot(aes(x = iter, y = nmi, color = method)) + 
#   # geom_line(aes(size = method), alpha = 0.5) +
#   geom_line(size = 1.2) + 
#   theme_minimal() +
#   # scale_colour_manual(values = c(1,1.5)) +
#   ggplot2::theme(
#     legend.background = ggplot2::element_blank(),
#     legend.title = ggplot2::element_blank(),
#     legend.position = c(0.8, 0.2),
#     # legend.text = ggplot2::element_text(size=18),
#   ) + 
#   ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
#   ylab("NMI") + xlab("Iteration") 

# hsbm::get_slice_nmi(xih_list2[[niter+1]], xi_tru)
# 
# xih_final1 = xih_list1[[niter+1]]
# xih_final2 = xih_list2[[niter+1]]
# 
# comp_all_matching_perms = function(xih) {
#   lapply(seq_along(xih), function(j) matching_perm(xih[[j]], xi_tru[[j]]))
# }
# comp_matching_score = function(xih_list){
#   length(unique(comp_all_matching_perms(xih_list))) # / length(xih_list)
# }
# 
# length(unique(comp_all_matching_perms(xih_final2)))
# length(unique(comp_all_matching_perms(xih_final1)))
# 
# comp_matching_score(xih_final2)
