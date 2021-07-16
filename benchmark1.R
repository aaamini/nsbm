Rcpp::sourceCpp("src/multsbm.cpp", verbose = T)
set.seed(1234)
n = 125
K = 7
alpha = 1
beta = 1
niter = 50
burnin = floor(niter/2)

z_tru = sample(1:K, n, replace = T)
B = nett::gen_rand_conn(n, K, lambda = 10)
A = nett::fast_sbm(z_tru, B)
zb = list(z_tru)

convert_cpp_label_matrix_to_list = function(zmat) {
  lapply(1:ncol(zmat), function(it) list(zmat[,it]+1))  # the list() is there to treat these  multi-layer labels with only a single layer 
}

Rcpp::sourceCpp("src/multsbm.cpp", verbose = T)

methods = list()
methods[["Mult-SBM-collaped-v1"]] = function(A) {
  z_list = convert_cpp_label_matrix_to_list(
    multsbm_collapsed_gibbs_sampler(A, K, alpha = alpha, beta = beta, niter = niter)
  )
  hsbm::get_map_labels(z_list, burnin = burnin, consecutive = T)$labels
}

methods[["Mult-SBM-collaped-v2"]] = function(A) {
  z_list = convert_cpp_label_matrix_to_list(
    multsbm_collapsed_gibbs_sampler_v2(A, K, alpha = alpha, beta = beta, niter = niter)
  )
  hsbm::get_map_labels(z_list, burnin = burnin, consecutive = T)$labels
}

methods[["Mult-SBM-regular"]] =  function(A) {
  z_list = convert_cpp_label_matrix_to_list(
    multsbm_gibbs_sampler_fast(A, K, alpha, beta = beta, niter = niter)
  )
  hsbm::get_map_labels(z_list, burnin = burnin, consecutive = T)$labels
}

mtd_names = names(methods)

res = do.call(rbind, lapply(seq_along(methods), function(j) {
  dt = as.numeric(system.time( zh <- methods[[j]](A) )["elapsed"])
  data.frame(method = mtd_names[j], 
             nmi = hsbm::get_agg_nmi(zb, zh), 
             elapsed_time = dt)
}))

print( knitr::kable(res, digits = 4, format="pipe") )