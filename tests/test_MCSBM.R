source("R/data_gen.R")
source("R/inference.R")
source("R/competing_methods.R")
source("R/matching.R")
source("R/plotting.R")
library(parallel)
library(ggplot2)
library(dplyr)
library(tibble)
library(Matrix)
library(RcppHungarian)
library(patchwork)

Rcpp::sourceCpp("src/MCSBM.cpp", verbose = T)
# Rcpp::sourceCpp("src/SBM.cpp", verbose = T)


# # Caompute the optimal matching permutation from zh to z, both label vectors
# matching_perm = function(zh, z) {
#   n = length(z)
#   out = HungarianSolver(1-compute_confusion_matrix(zh, z)/n)
#   perm = out$pairs[,2]
#   perm
# }
# 
# comp_all_matching_perms = function(xih, xi_tru) {
#   lapply(seq_along(xih), function(j) matching_perm(xih[[j]], xi_tru[[j]]))
# }
# comp_matching_score = function(xih, xi_tru){
#   length(unique(comp_all_matching_perms(xih, xi_tru))) # / length(xih)
# }

# set.seed(1235)
n = 100
K = 2
L = 3
J = 100

# lambda = NULL
lambda = 30
zeta = .9 # try z = 0.3
gam = .7
nreps = 10
n_cores = 32
niter = 300

nathan_data = F

if (nathan_data) { L = 3 }

methods = list()
# methods[["SBM"]] = function(A, L, niter=100) {
#   temp_label_list = lapply(seq_along(A), function(j) {
#     mod = new(SBM, A[[j]], L, 1, 1)
#     mod$run_gibbs(niter) 
#   })
#   lapply(1:(niter+1), function(iter) 
#     lapply(seq_along(temp_label_list), function(j) temp_label_list[[j]][,iter])
#   )
# }

methods[["MCSBM"]] = function(A, K, L, niter=50) {
#    out = mix_mcsbm(A, K, L, 1, 1, 0.9, 1, niter, 3)
  mod = new(MCSBM, A, K, L, 1, 1) #, 0.9, 0.99)
  out = mod$run_gibbs(niter)
#   
   list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
}

methods[["MCSBM-rnd"]] = function(A, K, L, niter=50) {
  #    out = mix_mcsbm(A, K, L, 1, 1, 0.9, 1, niter, 3)
  mod = new(MCSBM, A, K, L, 1, 1) #, 0.9, 0.99)
  mod$rnd_prob = 0.9
  mod$decay = .99
  out = mod$run_gibbs(niter)
  #   
  list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
}

methods[["MCSBM-mix"]] =  function(A, K, L, niter=50)  {
  out = mix_mcsbm(A, K, L, 1, 1, decay = 0.95, sa_temp = 100, rnd_prob = 0.5, niter, 3)  
  list(z = out$z, xi = lapply(out$xi, add_one_to_xi))
}

methods[["MCSBM-SA"]] = # fit_mcsbm
 function(A, K, L, niter=50) {
#  mod = new(MCSBM, A, K, L, 1, 1, 0, 0.95)
  mod = new(MCSBM, A, K, L, 1, 1)
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
    out = generate_nathans_data(n = n, J = J, K = K, lambda = lambda)  
  }
  else {
    out = gen_rand_nsbm(n=n, K=K, L=L, J=J, lambda=lambda, gam=gam, zeta=zeta, sort_z = T)
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
state_str =  sprintf("J = %d, n = %d, nr = %d, lam = %s, K = %d, L = %d", 
                     J, n, nreps, if(is.null(lambda)) "NA" else lambda, K, L)

if (nathan_data) {
  state_str = sprintf("%s_nathan", state_str)
} else {
  state_str = sprintf("%s, zeta = %1.1f, gam = %1.1f", state_str, zeta, gam)
}

p1 = plot_paths_and_avg(res, xi_nmi) + ylab("xi-NMI") + labs(title = state_str) 
p2 = plot_paths_and_avg(res, z_nmi, alpha_range = c(0.2, 1)) + ylab("z-NMI") # + labs(title = state_str)
p3 = plot_paths_and_avg(res, matching_score) + ylab("Matching score") # + labs(title = state_str)

print(p1 + p2 + p3)
# ggsave(sprintf("test_mcsbm2_%s.png", state_str), width = 10, height=5)

# p1 = res %>% 
#   group_by(iter, method) %>% summarise(xi_nmi = mean(xi_nmi)) %>% 
#   ggplot(aes(x = iter, y = xi_nmi, color = method)) + 
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
#   ylab("xi-NMI") + xlab("Iteration") + labs(title = state_str)
# 
# p2 = res %>% 
#   group_by(iter, method) %>% summarise(z_nmi = mean(z_nmi)) %>% 
#   ggplot(aes(x = iter, y = z_nmi, color = method)) + 
#   geom_line(size = 1.2) + 
#   # geom_violin() +
#   theme_minimal() +
#   # scale_colour_manual(values = c(1,1.5)) +
#   ggplot2::theme(
#     legend.background = ggplot2::element_blank(),
#     legend.title = ggplot2::element_blank(),
#     legend.position = c(0.8, 0.2),
#     # legend.text = ggplot2::element_text(size=18),
#   ) + 
#   ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
#   ylab("z-NMI") + xlab("Iteration") + labs(title = state_str)
# 
# p3 = res %>% 
#   group_by(iter, method) %>% summarise(matching_score = mean(matching_score)) %>% 
#   ggplot(aes(x = iter, y = matching_score, color = method)) + 
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
#   ylab("Matching score") + xlab("Iteration") + labs(title = state_str)
# 
# 
# p1 + p2 + p3
# # ggsave(sprintf("test_musbm_%s.png", state_str), width = 10, height=5)
# 
# # out = gen_rand_nsbm(n=n, K=K, L=L, J=J,  lambda=lambda, gam = gam, zeta=zeta, sort_z = T)
# # A = out$A
# # z_tru = out$z
# # xi_tru = out$xi
# # eta = out$eta[[1]]
# 
# # # check whether eta is permutation invariant
# # rand_perm = sample.int(L)
# # rand_P = as(rand_perm, "pMatrix")
# # sum((rand_P %*% eta %*% t(rand_P) - eta)^2)
# 
# 
# # # xih_list1 = lapply(seq_along(A), function(j) fit_sbm(A[[j]], L, niter))
# # xih_list1 = fit_sbm(A, L)
# # # xih_list1 = lapply(1:(niter+1), function(iter) 
# # #    lapply(seq_along(xih_list1), function(j) xih_list1[[j]][,iter])
# # #  )
# # 
# # 
# # xih_list2 = fit_musbm(A, L, niter)
# # # xih_list2 = lapply(xih_list2, add_one_to_xi)
# 
# # compute_nmi_evolution = function(xih_list) {
# #   sapply(seq_along(xih_list), function(iter) 
# #     hsbm::get_slice_nmi(xih_list[[iter]], xi_tru))
# # }
# # 
# # iter_seq = 1:(niter+1)
# # res = data.frame(iter = iter_seq, nmi = compute_nmi_evolution(xih_list1), method = "SBM")
# # res = rbind(res, 
# #             data.frame(iter = iter_seq, nmi = compute_nmi_evolution(xih_list2), method = "MuSBM")
# # )
# 
# # res %>% 
# #   ggplot(aes(x = iter, y = nmi, color = method)) + 
# #   # geom_line(aes(size = method), alpha = 0.5) +
# #   geom_line(size = 1.2) + 
# #   theme_minimal() +
# #   # scale_colour_manual(values = c(1,1.5)) +
# #   ggplot2::theme(
# #     legend.background = ggplot2::element_blank(),
# #     legend.title = ggplot2::element_blank(),
# #     legend.position = c(0.8, 0.2),
# #     # legend.text = ggplot2::element_text(size=18),
# #   ) + 
# #   ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
# #   ylab("NMI") + xlab("Iteration") 
# 
# # hsbm::get_slice_nmi(xih_list2[[niter+1]], xi_tru)
# # 
# # xih_final1 = xih_list1[[niter+1]]
# # xih_final2 = xih_list2[[niter+1]]
# # 
# # comp_all_matching_perms = function(xih) {
# #   lapply(seq_along(xih), function(j) matching_perm(xih[[j]], xi_tru[[j]]))
# # }
# # comp_matching_score = function(xih_list){
# #   length(unique(comp_all_matching_perms(xih_list))) # / length(xih_list)
# # }
# # 
# # length(unique(comp_all_matching_perms(xih_final2)))
# # length(unique(comp_all_matching_perms(xih_final1)))
# # 
# # comp_matching_score(xih_final2)
