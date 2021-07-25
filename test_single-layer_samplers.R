Rcpp::sourceCpp("src/multsbm.cpp", verbose = T)

library(ggplot2)
# theme_set(theme_minimal(base_size = 18))
theme_set(theme_minimal())
library(dplyr)

set.seed(185)
n = 50 # try 125
lambda = 10 # SBM average expected degree 
K = 3
alpha = 1
beta = 1
Zcap = 5  # increasing this slows down the collapsed one significantly

niter = 50
nreps = 50
nreps_per_net = 3 # Try also 5
include_dpsbm_flag = F

comp_agg_nmi_path = function(z_list) {
   sapply(seq_along(z_list), function(it) hsbm::get_agg_nmi(z_list[[it]], list(z_tru))) 
}

convert_cpp_label_matrix_to_list = function(zmat) {
  lapply(1:ncol(zmat), function(it) list(zmat[,it]+1))  # the list() is there to treat these  multi-layer labels with only a single layer 
}

methods = list()
# methods[["Mult-SBM-collaped-v1"]] = function(A) { 
#   multsbm_collapsed_gibbs_sampler(A, K, alpha = alpha, beta = beta, niter = niter)
# }

methods[["Mult-SBM-collaped-v2"]] = function(A) {
  multsbm_collapsed_gibbs_sampler_v2(A, K, alpha = alpha, beta = beta, niter = niter)
}

methods[["Mult-SBM-collaped-v3"]] = function(A) {
  multsbm_collapsed_gibbs_sampler_v3(A, K, alpha = alpha, beta = beta, niter = niter)
}

methods[["Mult-SBM-collaped-v4"]] = function(A) {
  multsbm_collapsed_gibbs_sampler_v4(A, K, alpha = alpha, beta = beta, niter = niter)
}


methods[["Mult-SBM-regular"]] =  function(A) {
  multsbm_gibbs_sampler_fast(A, K, alpha, beta = beta, niter = niter)
}

if (include_dpsbm_flag) {
  # DP-SBM
  Rcpp::sourceCpp("src/dpsbm.cpp", verbose = T)

  methods[["DP-SBM-collapsed-v1"]] =  function(A) {
    fit_dpsbm_collapsed(A, Zcap = Zcap, niter = niter)
  }

  methods[["DP-SBM-regular"]] =  function(A) {
    fit_dpsbm(A, Zcap = Zcap, niter = niter)
  }
}


mtd_names = names(methods)

res = NULL
for (rep in 1:nreps) {
  if (!rep %% 5) cat('.')
  z_tru = sample(1:K, n, replace = T)
  B = nett::gen_rand_conn(n, K, lambda = lambda) #, gamma = 0.5)
  A = nett::fast_sbm(z_tru, B)
  # diag(A) = sample(0:1, n, replace = T)  
  
  for (j in 1:nreps_per_net) {
    res_curr = do.call(rbind, lapply(seq_along(methods), function(j) {
       z_list = convert_cpp_label_matrix_to_list( methods[[j]](A) )
       data.frame(iter = 1:niter,
          rep = rep,
          rep_per_net = j,
          nmi = comp_agg_nmi_path(z_list),
          method = mtd_names[j])
    }))

    res = rbind(res, res_curr)
  }
}



res %>% 
  group_by(iter, method) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = iter, y = nmi, color = method)) + 
  geom_line() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("Average NMI") + xlab("Iteration") +
  labs(title = sprintf("n = %d,  lambda = %2.1f", n, lambda))
ggsave(sprintf("collapsed_vs_regular_n=%d_lam=%d.png", n, round(lambda)), 
    width = 5, height = 4)

# Tests 
# hsbm::seq_nmi_plot(z_list)
# comp_beta_matrix(A, z-1, K, alpha1, beta1)
# comp_Bet(A, z, K, alpha1, beta1)

# # Speed test
# microbenchmark::microbenchmark(regular =  fit_dpsbm(A, Zcap = Zcap, niter = niter),
#                                 collapsed =  fit_dpsbm_collapsed(A, Zcap = Zcap, niter = niter), times = 20)
