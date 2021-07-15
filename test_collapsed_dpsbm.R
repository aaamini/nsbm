Rcpp::sourceCpp("src/dpsbm.cpp", verbose = T)
set.seed(185)
n = 50
K = 3
# alpha1 = 1
# beta1 = 1
Zcap = 5  # increasing this slows down the collapsed one significantly

niter = 50
nreps = 50
nreps_per_net = 3 # Try also 5

comp_agg_nmi_path = function(z_list) {
   sapply(seq_along(z_list), function(it) hsbm::get_agg_nmi(z_list[[it]], list(z_tru))) 
}

convert_cpp_label_matrix_to_list = function(zmat) {
  lapply(1:ncol(zmat), function(it) list(zmat[,it]+1))  # the list() is there to treat these  multi-layer labels with only a single layer 
}

res = NULL
for (rep in 1:nreps) {
  z_tru = sample(1:K, n, replace = T)
  B = nett::gen_rand_conn(n, K, lambda = 10)
  A = nett::fast_sbm(z_tru, B)
  # diag(A) = sample(0:1, n, replace = T)  
  
  for (j in 1:nreps_per_net) {
    # C++ version of the collapsed sampler
    z_list = convert_cpp_label_matrix_to_list(
      fit_dpsbm_collapsed(A, Zcap = Zcap, niter = niter)
    )
    res = rbind(res, data.frame(iter = 1:niter,
                                rep = rep,
                                rep_per_net = j,
                                nmi = comp_agg_nmi_path(z_list),
                                method = "collapsed"))
    
    z_list = convert_cpp_label_matrix_to_list(
      fit_dpsbm(A, Zcap = Zcap, niter = niter)
    )
    # z_list = lapply(1:niter, function(it) list(out[,it]+1))
    res = rbind(res, data.frame(iter = 1:niter, 
                                rep = rep, 
                                rep_per_net = j, 
                                nmi = comp_agg_nmi_path(z_list), 
                                method = "regular"))
  }
}


library(ggplot2)
library(dplyr)
res %>% 
  group_by(iter, method) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = iter, y = nmi, color = method)) + 
  geom_line() + theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.1),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = 1.25)) +
  ylab("Average NMI") + xlab("Iteration")
ggsave("collapsed_vs_regular.png", width = 5, height = 4)

# Tests 
# hsbm::seq_nmi_plot(z_list)
# comp_beta_matrix(A, z-1, K, alpha1, beta1)
# comp_Bet(A, z, K, alpha1, beta1)

# # Speed test
microbenchmark::microbenchmark(regular =  fit_dpsbm(A, Zcap = Zcap, niter = niter),
                                collapsed =  fit_dpsbm_collapsed(A, Zcap = Zcap, niter = niter), times = 20)
