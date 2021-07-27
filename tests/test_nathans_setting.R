# library ---
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())
source("R/nSBM_functions.R")
source("R/splice_sampler.R")
source("R/data_gen.R")

# simulation ----
# set.seed(575)
niter = 10
K = L = 10
n_cores = 32
sparse_data = F

if (!sparse_data) {
  n = 40; J= 60
} else{
  n = 50; J = 25; lambda = 10; K_tru = 3
}

nreps = 12
# res = do.call(rbind, lapply(1:nreps, function(rep) {
res = do.call(rbind, mclapply(1:nreps, function(rep) {
  if (!sparse_data){
    out = generate_nathans_data(n = n, J = J) 
  } else {
    out = generate_sparse_random_data(n = n, J = J, lambda = lambda, K_tru = K_tru) 
  }

  A = out$A
  z_tru = out$z
  
  res = NULL
  
  dt = system.time(samp <- gibbs.nSBM(A, K, L, ns = niter, monitor = TRUE))["elapsed"]
  res = rbind(res, data.frame(
    rep = rep,
    dt = dt,
    iter = 1:(niter + 1), 
    nmi = apply(samp$z, 2, function(z) nett::compute_mutual_info(z, z_tru)),
    method = "Nathan's"))
  
  dt = system.time(samp <- splice_sampler(A, K, L, ns = niter, monitor = TRUE))["elapsed"]
  res = rbind(res, data.frame(
    rep = rep,
    dt = dt,
    iter = 1:(niter + 1), 
    nmi = apply(samp$z, 2, function(z) nett::compute_mutual_info(z, z_tru)),
    method = "Splice"))
  
  dt = system.time({
    model =  new(NestedSBM, A, K, L)
    fitted_model <- model$run_gibbs(niter) 
  })["elapsed"]
  res = rbind(res, data.frame(
    rep = rep,
    dt = dt,
    iter = 1:(niter + 1), 
    nmi = apply(fitted_model$z, 2, function(z) nett::compute_mutual_info(z, z_tru)),
    method = "Pure C++"))
# }))
}, mc.cores = n_cores))


state_str =  sprintf("J = %d, n = %d, nreps = %d", J, n, nreps)
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
  ylab("NMI") + xlab("Iteration") + labs(title = state_str)


tag = "sparse"
ggsave(sprintf("test_splice_new_%s_%s.png", state_str, tag), width = 6, height=5)

knitr::kable( res %>% 
  filter(iter > niter / 2) %>% 
  group_by(method) %>% summarise(dt = mean(dt), nmi = mean(nmi)) 
)
