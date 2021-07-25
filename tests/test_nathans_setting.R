# library ---
library(igraph)
library(ggplot2)
library(dplyr)
Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

# simulation ----
set.seed(575)
m <- 60
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
niter = 10
K = L = 10
model = new(NestedSBM, A, K, L)
z_init = model$z + 1
system.time( fitted_model <- model$run_gibbs(niter) )["elapsed"]

res = data.frame(iter = 1:(niter + 1), 
                 nmi = apply(cbind(z_init, fitted_model$z), 2, function(z) nett::compute_mutual_info(z, z_tru)),
                 method = "NSBM (C++)")

source("nathans_nsbm/nSBM_functions.R")
system.time(samp <- gibbs.nSBM(A, K, L, ns = niter, monitor = TRUE))["elapsed"]

# z_init is not relevant to this, but added to make the sequence of the same size.
res = rbind(res, data.frame(
            iter = 1:(niter + 1), 
            nmi = apply(cbind(z_init, samp$z), 2, function(z) nett::compute_mutual_info(z, z_tru)),
            method = "NSBM (Nathan's)"))

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

ggsave("test.png", width = 6, height=5)

