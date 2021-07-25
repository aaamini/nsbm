# library ---
library(igraph)
library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

# simulation ----
# set.seed(575)
m <- 6
pm1 <- cbind(c(.9, .75, .5)
             , c(.75, .6, .25)
             , c(.5, .25, .1))
pm2 <- cbind(c(.8, .1, .3)
             , c(.1, .9, .2)
             , c(.3, .2, .7))
pm3 <- cbind(c(.1, .4, .6)
             , c(.4, .3, .1)
             , c(.6, .1, .5))

K.true <- c(rep(1, m/3), rep(2, m/3), rep(3, m/3))
A <- L.true <- vector("list", length = m)
for (j in seq_len(m)) {
  
  # vary network order
  # n <- sample(seq(60, 80, 20), 1)
  n <- 80
  
  if (j <= m/3) {
    G <- as_adj(sample_sbm(n, pm1, c(.4*n, .35*n, .25*n)))
    L.true[[j]] <- c(rep(1, .4*n), rep(2, .35*n), rep(3, .25*n))
  
  } else if (j <= 2*m/3) {
    G <- as_adj(sample_sbm(n, pm2, c(.7*n, .15*n, .15*n)))
    L.true[[j]] <- c(rep(1, .7*n), rep(2, .15*n), rep(3, .15*n))
    
  } else {
    G <- as_adj(sample_sbm(n, pm3, c(.2*n, .4*n, .4*n)))
    L.true[[j]] <- c(rep(1, .2*n), rep(2, .4*n), rep(3, .4*n))
  }
  
  # unlabeled networks
  pi <- sample(n)
  A[[j]] <- G[pi, pi]
  L.true[[j]] <- L.true[[j]][pi]
}

z_tru = K.true
xi_tru = L.true
niter = 20
K = L = 10
model = new(NestedSBM, A, K, L)
z_init = model$z + 1
xi_init = lapply(model$xi, function(x) x + 1)
system.time( fitted_model <- model$run_gibbs(niter) )["elapsed"]
# system.time( fitted_model <- model$run_gibbs_naive(niter) )["elapsed"]

xi_hist = lapply(fitted_model$xi, function(xi) lapply(xi, function(x) x + 1))

compute_xi_nmi = function(xi_hist) {
  sapply(c(list(xi_init), xi_hist), function(xi) hsbm::get_slice_nmi(xi, xi_tru))
}
mat_to_list = function(x) lapply(seq_len(ncol(x)), function(i) x[,i])

res = data.frame(iter = 1:(niter + 1), 
                 nmi = compute_xi_nmi(xi_hist),
                 method = "NSBM (C++)")

source("nathans_nsbm/nSBM_functions.R")
system.time(samp <- gibbs.nSBM(A, K, L, ns = niter, monitor = TRUE))["elapsed"]

# z_init is not relevant to this, but added to make the sequence of the same size.
res = rbind(res, data.frame(
            iter = 1:(niter + 1), 
            nmi = compute_xi_nmi(lapply(1:niter, function(iter) mat_to_list(samp$xi[,,iter]))),
            method = "NSBM (Nathan's)"))

model$z
nett::compute_mutual_info(model$xi[[3]]+1, xi_tru[[3]])
cbind(model$xi[[3]]+1, xi_tru[[3]])

nett::compute_mutual_info(as.vector(model$z)+1, z_tru)

res %>% 
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
  ylab("NMI") + xlab("Iteration") 

ggsave("test3.png", width = 6, height=5)


