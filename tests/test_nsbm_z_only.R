# Clean test file -- should run without trouble

Rcpp::sourceCpp("src/utils.cpp", verbose = T)
Rcpp::sourceCpp("src/nsbm.cpp", verbose = T)
#set.seed(1234)
n = 30
K = 4
L = 5
J = 25
alpha = 1
beta = 1
niter = 5
lambda = 7

# data
z_tru = sample(1:K, J, replace = T)
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = 7))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

# tensor comp
comp_m_tensor = function(A, z, xi, K, L) {
    lapply(1:K, function(k) {
      idx = which(z == k)
      if (length(idx) == 0) return( matrix(0, L, L) )
      Reduce(`+`, lapply(idx, 
                         function(j) comp_blk_sums_and_sizes(A[[j]], xi[[j]]-1, L)$lambda))  
    }
  )
}

comp_mbar_tensor = function(A, z, xi, K, L) {
  lapply(1:K, function(k) {
    idx = which(z == k)
    if (length(idx) == 0) return( matrix(0, L, L) )
    Reduce(`+`, lapply(idx,  
                       function(j) {
                         out = comp_blk_sums_and_sizes(A[[j]], xi[[j]]-1, L)
                         out$NN - out$lambda
                       }))
    }
  )
}

upper_tri_sum = function(X) {
  sum(X - diag(diag(X)))/2 + sum(diag(X))
}

upper_tri_sum2 = function(X) {
  uidx = upper.tri(diag(nrow(X)), diag = T)
  sum(X[uidx])
}

remove_max = function(x) x - max(x)

## tests
z_init = sample(1:K, J, replace = T) # z_tru
# m_tensor = comp_m_tensor(A, z_prev, xi_tru, K,  L)
# mbar_tensor = comp_mbar_tensor(A, z_prev, xi_tru, K,  L)
# q = simplify2array(m_tensor) # turn list into 3d array
# qbar = simplify2array(mbar_tensor)

# out = comp_blk_sums_and_sizes(A[[j]], xi_tru[[j]]-1, L)
# D = out$lambda
# Dbar = out$NN - D

nsbm_z_update_R = function(A, z_init, xi_tru, K, L, alpha, beta, niter = 5) {
  z_prev = z_init
  z_list = list(z_prev)
  for (iter in 1:niter) {
    for (j in 1:J) {
      # out = nsbm_single_z_update_R(A, z_prev, j, xi_tru, K, L, alpha, beta)
      log_prob_R = rep(0, K)
      r0 = z_prev[j]
      z = z_prev
      for (r in 1:K) {
        z[j] = r
        q_new = simplify2array(comp_m_tensor(A, z, xi_tru, K,  L))
        qbar_new = simplify2array(comp_mbar_tensor(A, z, xi_tru, K,  L))
        
        #temp = log(beta(q_new + alpha, qbar_new + beta)) - 
        #  log(beta(q + alpha, qbar + beta))
        temp = log(beta(q_new + alpha, qbar_new + beta))
        temp_k_sum = rowSums(temp, dims=2)
        log_prob_R[r] = upper_tri_sum(temp_k_sum) 
      }
      z[j] = sample(K, 1, prob = safe_exp(log_prob_R))
      # print(out$log_prob)
      z_list = c(z_list, list(z))
      z_prev = z
    }
  }
  z_list
  # list(z=z, log_prob = remove_max(log_prob_R))
}

# out = nsbm_single_z_update(A, z_prev, j, xi_tru, K, L, alpha, beta)
z_list = nsbm_z_update_R(A, z_init, xi_tru, K, L, alpha, beta, niter = niter)

res = data.frame(inner_iter = 1:(niter*J + 1), 
                 nmi = sapply(z_list, function(z) nett::compute_mutual_info(z, z_tru)),
                 method = "R")



xi_tru_0idx = lapply(xi_tru, function(xi) xi-1)
# m = mbar = array(0, dim = c(L,L,K))
# comp_count_tensors(A, z_tru-1, xi_tru_0idx, L, m, mbar)

z_mat_cpp = nsbm_z_update_cpp(A, xi_tru_0idx, z_init-1, L = L, K = K, niter = 5)
res =  rbind(res, 
             data.frame(inner_iter = 1:(niter*J + 1), 
                        nmi = apply(cbind(z_init, z_mat_cpp), 2, function(z) nett::compute_mutual_info(z, z_tru)),
                        method = "C++"))

library(ggplot2)
library(dplyr)
res %>% 
  ggplot(aes(x = inner_iter, y = nmi, color = method)) + 
  geom_line(aes(size = method), alpha = 0.5) + theme_minimal() +
  # scale_colour_manual(values = c(1,1.5)) +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("NMI") + xlab("Inner Iteration") 

ggsave(sprintf("nsbm_z_update_test_n=%d_J=%d.png", n, round(lambda)), 
       width = 5, height = 4)
# q = out$m
# qbar = out$mbar
# q_R = simplify2array(comp_m_tensor(A, z_tru, xi_tru, K,  L))
# qbar_R = simplify2array(comp_mbar_tensor(A, z_tru, xi_tru, K,  L))
# 
# sum(abs(q-q_R))
# sum(abs(qbar - qbar_R))


# Rcpp::sourceCpp("src/utils.cpp", verbose = T)
# 
# log_prob_cpp = as.vector(
#   comp_tensor_log_beta_ratio_sums(q, qbar, D, Dbar, r0-1, alpha, beta)
# )

#remove_max = function(x) x - max(x)
#rbind(remove_max(log_prob_R), remove_max(log_prob_cpp))


# rbind(log_prob_R, log_prob_cpp)
# cat(sprintf('Error = %e', sum(abs(log_prob_R - log_prob_cpp))))
# 
# c(r0, which.max(log_prob_cpp))
