# library ---
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
library(patchwork)
#Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
#setMethod("show", "Rcpp_NestedSBM", function(object) object$print())
source("R/nSBM_functions.R")
source("R/splice_sampler.R")
source("R/data_gen.R")
source("R/inference.R")
source("R/nsbm_wrapper.R")

# simulation ----
set.seed(1234)
niter = 300
K = L = 10
n_cores = 12
nreps = 10
sparse_data = T

if (!sparse_data) {
  n = 50; J = 12
} else{
  n = 100; J = 16; lambda = 15; K_tru = 3
}


methods = list()
methods[["C++ (collapsed)"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = T)
  # model =  new(NestedSBM, A, K, L)
  # model$run_gibbs(niter)
}

methods[["Splice"]] = function(A) {
  splice_sampler(A, K, L, ns = niter, monitor = TRUE)
}

methods[["C++ (non-collapsed v1)"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 1)
  # model =  new(NestedSBM, A, K, L)
  # model$run_gibbs_via_eta(niter) 
}

methods[["C++ (non-collapsed v2)"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 2)
  # model =  new(NestedSBM, A, K, L)
  # model$run_gibbs_via_eta(niter) 
}

methods[["C++ (non-collapsed v3)"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 3)
  # model =  new(NestedSBM, A, K, L)
  # model$run_gibbs_via_eta(niter) 
}


# methods[["Nathan's"]] = function(A) {
#   samp <- gibbs.nSBM(A, K, L, ns = niter, monitor = TRUE)
#   samp$z
# }

mtd_names = names(methods)

res = do.call(rbind, mclapply(1:nreps, function(rep) {
# res = do.call(rbind, lapply(1:nreps, function(rep) {
    if (!sparse_data){
        out = generate_nathans_data(n = n, J = J) 
    } else {
        # out = generate_sparse_random_data(n = n, J = J, lambda = lambda, K_tru = K_tru) 
      out = gen_rand_nsbm(n=n, K = K_tru, J = J, gam = 0.3, zeta=0.5, lambda = lambda)
    }
    A = out$A
    z_tru = out$z
    xi_tru = out$xi

    do.call(rbind, lapply(seq_along(methods), function(j) { 
        dt = as.numeric(system.time( mout <- methods[[j]](A) )["elapsed"])
        z_hist = mout$z
        xi_hist = mout$xi
        data.frame(method = mtd_names[j], 
                    rep = rep,
                    iter = 1:(niter + 1), 
                    z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru)),
                    xi_nmi = sapply(xi_hist, function(xi) hsbm::get_slice_nmi(xi, xi_tru)),
                    elapsed_time = dt)
    }))
}, mc.cores = n_cores))
# }))    


state_str =  sprintf("J = %d, n = %d, nreps = %d", J, n, nreps)
p1 = res %>% 
  group_by(iter, method) %>% summarise(z_nmi = mean(z_nmi)) %>% 
  ggplot(aes(x = iter, y = z_nmi, color = method)) + 
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
  ylab("z-NMI") + xlab("Iteration") + labs(title = state_str)

# print(p1)

# tag = "sparse"
#tag = "nathan"
#ggsave(sprintf("test_non_collapsed_z_nmi_%s_%s.png", state_str, tag), width = 6, height=5)

p2 = res %>% 
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

# print(p2)

print(p1 + p2)
tag = "sparse"
#tag = "nathan"
ggsave(sprintf("test_many_v_%s_%s.png", state_str, tag), width = 10, height=5)

#ggsave(sprintf("test_non_collapsed_xi_nmi_%s_%s.png", state_str, tag), width = 6, height=5)

print(knitr::kable( res %>% 
                      filter(iter > niter / 2) %>% 
                      group_by(method) %>% summarise(elapsed_time = mean(elapsed_time), 
                                                     z_nmi = mean(z_nmi),
                                                     xi_nmi = mean(xi_nmi)) 
))
