library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/models.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

# set.seed(1234)
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

## Set z to true z and update xi only
model$z = z_tru-1
model$z 
model$set_xi_to_random_labels()
xi_init = lapply(model$xi, function(x) x + 1)
xi_init
model$comp_count_tensors()
model$m

niter = 200
xi_list = vector("list", niter*n)
j = 5
for (iter in 1:niter - 1) {
  model$update_w()
  # model$update_pi()
  for (j in 1:J - 1) {
    for (s in 1:n - 1) {
      model$update_xi_element(j, s)
      xi_list[[n*iter + s + 1]] = as.vector(model$xi[[j+1]]) + 1
    }
    # model$update_z_element(j)
  }
}

# nett::compute_mutual_info(as.vector(model$z)+1, z_tru)

res = NULL
res =  rbind(res, 
             data.frame(inner_iter = 1:(niter*n + 1), 
                        nmi = apply(cbind(xi_init[[j]], simplify2array(xi_list)), 2, function(xi) nett::compute_mutual_info(xi, xi_tru[[j]])),
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



# model$xi = xi_tru
# model$xi[[10]] - xi_tru[[10]]
