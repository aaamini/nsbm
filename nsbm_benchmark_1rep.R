# library ---
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
library(patchwork)

# setwd("/project/sand/njosephs/NDP/nsbm/")
#Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
#setMethod("show", "Rcpp_NestedSBM", function(object) object$print())
source("R/nSBM_functions.R")
source("R/splice_sampler.R")
source("R/data_gen.R")
source("R/inference.R")
source("R/nsbm_wrapper.R")

# simulation ----
#set.seed(1234)
niter = 100
K = L = 10
n_cores <- 1# detectCores()
nreps <- n_cores
rand_sbm <- TRUE
tag <- "A"

n = 100; J = 20; K_tru = 3; L_tru = c(2,3,5)

methods = list()

methods[["Collapsed"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = T, version = 1, naive = TRUE)
}

methods[["Gibbs"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 3)
}

methods[["Block Gibbs"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 5)
}

methods[["Incompatible block"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 6)
}

# methods[["Block collapsed"]] = function(A) {
#   fit_nsbm(A, K, L, niter, collapsed = T, version = 2, naive = TRUE)
# }

mtd_names = names(methods)

res = do.call(rbind, lapply(1:nreps, function(rep) {
  
  if (rand_sbm) {
    out = gen_rand_nsbm(n = n
                        , J = J
                        , K = K_tru
                        , L = L_tru
                        , gam = 0.4, lambda = 15)    
  } else {
    out = generate_nathans_data()  
  }
  
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
  
  do.call(rbind, lapply(seq_along(methods), function(j) {
    cat("Method:", mtd_names[j], "...")
    dt = as.numeric(system.time( mout <- methods[[j]](A) )["elapsed"])
    z_hist = mout$z
    xi_hist = mout$xi
    
    cat(sprintf("%3.2f\n", dt))
    data.frame(method = mtd_names[j], 
               rep = rep,
               iter = 1:(niter + 1), 
               z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru)),
               xi_nmi = sapply(xi_hist, function(xi) hsbm::get_slice_nmi(xi, xi_tru)),
               elapsed_time = dt) 
  }))
}))

res$method <- factor(res$method, levels = c("Gibbs"
                                            , "Collapsed"
                                            , "Block Gibbs"
                                            , "Incompatible block"
                                            , "Block collapsed"))

if (rand_sbm) {
  state_str =  sprintf("J = %d, n = %d, nreps = %d", J, n, nreps)  
} else {
  state_str = "HSBM: Multilayer personality-friendship network"
}

p1 = res %>% 
  group_by(iter, method) %>% summarise(mean = mean(z_nmi), sd = sd(z_nmi)) %>% 
  ggplot(aes(x = iter, y = mean, color = method)) + 
  geom_line(size = 1.4) + 
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("z-NMI") + xlab("Iteration") + labs(title = state_str)

p2 = res %>% 
  group_by(iter, method) %>% summarise(mean = mean(xi_nmi), sd = sd(xi_nmi)) %>% 
  ggplot(aes(x = iter, y = mean, color = method)) + 
  geom_line(size = 1.4) + 
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("xi-NMI") + xlab("Iteration")

print(p1 + p2)

ggsave(sprintf("test_many_v_%s_%s.png", state_str, tag), width = 10, height=5)

print(knitr::kable( res %>% 
                      filter(iter > niter / 2) %>% 
                      group_by(method) %>% summarise(elapsed_time = mean(elapsed_time), 
                                                     z_nmi = mean(z_nmi),
                                                     xi_nmi = mean(xi_nmi)) 
))
