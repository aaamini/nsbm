library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())


#set.seed(1234)
n = 30
K = 4
L = 5
J = 20
alpha = 1
beta = 1
niter = 10
lambda = 7

# data
z_tru = sample(1:K, J, replace = T)
eta = lapply(1:K, function(i) nett::gen_rand_conn(n, L, lambda = lambda))
xi_tru = A = vector("list", J)
for (j in 1:J) {
  xi_tru[[j]] = sample(L, n, replace = T, )
  A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
}

z_init = sample(1:K, J, replace = T) #
z_init_cpp = z_init - 1
xi_tru_cpp = lapply(xi_tru, function(x) x-1)


model = new(NestedSBM, A, K, L)
model

## Set xi to true xi and update z only
model$z = z_init_cpp
model$xi = xi_tru_cpp
model$comp_count_tensors()

z_list = vector("list", niter*J)
for (iter in 1:niter - 1) {
  for (j in 1:J - 1) {
    model$update_z_element(j)
    z_list[[J*iter + j + 1]] = as.vector(model$z) + 1
  }
}

res = NULL
res =  rbind(res, 
             data.frame(inner_iter = 1:(niter*J + 1), 
                        nmi = apply(cbind(z_init, simplify2array(z_list)), 2, function(z) nett::compute_mutual_info(z, z_tru)),
                        method = "NSBM (C++)"))
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
