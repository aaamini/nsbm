# library ---
library(igraph)
library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

# simulation ----
set.seed(575)
niter = 10
K = L = 10
n = 40
m = 6

source("nathans_nsbm/data_gen.R")
# model = new(NestedSBM, A, K, L)
# z_init = model$z + 1
# # system.time( fitted_model <- model$run_gibbs(niter) )["elapsed"]
# system.time( fitted_model <- model$run_gibbs_naive(niter) )["elapsed"]
# 
# res = data.frame(iter = 1:(niter + 1), 
#                  nmi = apply(cbind(z_init, fitted_model$z), 2, function(z) nett::compute_mutual_info(z, z_tru)),
#                  method = "NSBM (C++)")

source("nathans_nsbm/nSBM_functions.R")
source("nathans_nsbm/spliced_sampler.R")

nreps = 2
res = do.call(rbind, lapply(1:nreps, function(rep) {
  out = generate_nathans_data(n = n, m = m)
  A = out$A
  z_tru = out$z
  
  res = NULL
  system.time(samp <- gibbs.nSBM(A, K, L, ns = niter, monitor = TRUE))["elapsed"]
  # z_init is not relevant to this, but added to make the sequence of the same size.
  res = rbind(res, data.frame(
    rep = rep,
    iter = 1:(niter + 1), 
    nmi = apply(samp$z, 2, function(z) nett::compute_mutual_info(z, z_tru)),
    method = "NSBM (Nathan's)"))
  
  system.time(samp <- splice_sampler(A, K, L, ns = niter, monitor = TRUE))["elapsed"]
  
  res = rbind(res, data.frame(
    rep = rep,
    iter = 1:(niter + 1), 
    nmi = apply(samp$z, 2, function(z) nett::compute_mutual_info(z, z_tru)),
    method = "Splice"))
}))

res %>% 
  group_by(iter, method) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = iter, y = nmi, color = method)) + 
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
  ylab("NMI") + xlab("Iteration") + labs(title = sprintf("m = %d, n = %d", m, n))

# ggsave("test_splice.png", width = 6, height=5)

