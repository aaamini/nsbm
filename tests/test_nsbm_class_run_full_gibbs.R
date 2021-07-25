library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

# set.seed(1234)
n = 30
K = 4
L = 5
J = 20
alpha = 1
beta = 1
niter = 100
lambda = 7

# data
z_tru = sample(1:K, J, replace = T)
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = lambda))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

model = new(NestedSBM, A, K, L)
z_init = model$z + 1
fitted_model = model$run_gibbs(niter)


res = NULL
res =  rbind(res, 
             data.frame(inner_iter = 1:(niter + 1), 
                        nmi = apply(cbind(z_init, fitted_model$z), 2, function(z) nett::compute_mutual_info(z, z_tru)),
                        method = "NSBM (C++)"))

hsbm::seq_nmi_plot(lapply(1:niter, function(iter) list(fitted_model$z[,iter])))

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
