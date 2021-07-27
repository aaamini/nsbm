library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

# set.seed(1234)
n = 50
K = 10
L = 10
J = 25
alpha = 1
beta = 1
niter = 100
lambda = 10

# data
K_tru = 3


res = NULL

# model = new(NestedSBM, A, K, L)
# z_init = model$z + 1
# fitted_model = model$run_gibbs(niter)
# res =  rbind(res, 
#              data.frame(inner_iter = 1:(niter + 1), 
#                         nmi = apply(cbind(z_init, fitted_model$z), 2, function(z) nett::compute_mutual_info(z, z_tru)),
#                         method = "NSBM (C++)"))

nreps = 1
res = do.call(rbind, lapply(1:nreps, function(rep) {
  # res = do.call(rbind, mclapply(1:nreps, function(rep) {
  z_tru = sample(1:K_tru, J, replace = T)
  eta = lapply(1:K_tru, function(i) nett::gen_rand_conn(n, L, lambda = lambda))
  xi_tru = A = vector("list", J)
  for (j in 1:J) {
    xi_tru[[j]] = sample(L, n, replace = T, )
    A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
  }
  
  res = NULL
  
  dt = system.time(samp <- gibbs.nSBM(A, K, L, ns = niter, monitor = TRUE))["elapsed"]
  res = rbind(res, data.frame(
    rep = rep,
    dt = dt,
    iter = 1:(niter + 1), 
    nmi = apply(samp$z, 2, function(z) nett::compute_mutual_info(z, z_tru)),
    method = "NSBM (Nathan's)"))
  
  dt = system.time(samp <- splice_sampler(A, K, L, ns = niter, monitor = TRUE))["elapsed"]
  res = rbind(res, data.frame(
    rep = rep,
    dt = dt,
    iter = 1:(niter + 1), 
    nmi = apply(samp$z, 2, function(z) nett::compute_mutual_info(z, z_tru)),
    method = "Splice"))
}))
hsbm::seq_nmi_plot(lapply(1:niter, function(iter) list(fitted_model$z[,iter])))

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
  ylab("NMI") + xlab("Iteration") + labs(title = sprintf("m = %d, n = %d, nreps = %d", m, n, nreps))

# ggsave("test_splice3.png", width = 6, height=5)
# ggsave("test_splice_arash_setting.png", width = 6, height=5)

res %>% 
  group_by(method) %>% summarise(dt = mean(dt)) 

# # R version of model$run_gibbs()
# niter = 100
# z_hist = xi_hist = vector("list", niter)
# model$comp_count_tensors()  # init count tensors
# for (iter in 1:niter - 1) {
#   model$update_w()
#   model$update_pi()
#   for (j in 1:J - 1) {
#     for (s in 1:n - 1) {
#       model$update_xi_element(j, s)
#     }
#     model$update_z_element(j)
#   }
#   xi_hist[[iter+1]] = model$xi
#   z_hist[[iter+1]] = model$z
# }
