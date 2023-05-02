# Libraries ----
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
library(patchwork)

# Functions ----
source("R/data_gen.R")
source("R/inference.R")
source("R/nsbm_wrapper.R")
source("R/NCLM.R")
source("R/setup_methods.R")

sample_hsbm <- function(J, n, labeled = TRUE) {
  
  pm1 <- cbind(c(.9, .75, .5)
               , c(.75, .6, .25)
               , c(.5, .25, .1))
  pm2 <- cbind(c(.8, .1, .3)
               , c(.1, .9, .2)
               , c(.3, .2, .7))
  pm3 <- cbind(c(.1, .4, .6)
               , c(.4, .3, .1)
               , c(.6, .1, .5))
  
  z_tru <- c(rep(1,J/3), rep(2,J/3), rep(3,J/3))
  A <- xi_tru <- vector("list", length = J)
  for (j in seq_len(J)) {
    
    if (j <=J/3) {
      G <- as_adj(sample_sbm(n, pm1, c(.4*n, .35*n, .25*n)))
      xi_tru[[j]] <- c(rep(1, .4*n), rep(2, .35*n), rep(3, .25*n))
      
    } else if (j <= 2*J/3) {
      G <- as_adj(sample_sbm(n, pm2, c(.7*n, .15*n, .15*n)))
      xi_tru[[j]] <- c(rep(1, .7*n), rep(2, .15*n), rep(3, .15*n))
      
    } else {
      G <- as_adj(sample_sbm(n, pm3, c(.2*n, .4*n, .4*n)))
      xi_tru[[j]] <- c(rep(1, .2*n), rep(2, .4*n), rep(3, .4*n))
    }
    
    if (labeled) {
      A[[j]] <- G
    } else {
      pi <- sample(n)
      A[[j]] <- G[pi, pi]
      xi_tru[[j]] <- xi_tru[[j]][pi]
    }
    
  }
  
  return(list(A = A, z = z_tru, xi = xi_tru))
}

# Settings ----
niter <- 200  # number of iteration for Gibbs samplers
K <- L <- 15  # truncation levels for NSBM models
ncores <- detectCores()
nreps <- 100

J <- 120 # number of networks
n <- 60  # number of nodes

# Simulation ----
res <- do.call(rbind, mclapply(seq_len(nreps), function(rep) {
  set.seed(rep)
  
  out <- sample_hsbm(J, n)
  
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
  
  out = do.call(rbind, lapply(seq_along(methods), function(j) {
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    end_time = as.numeric(Sys.time() - start_time)
    
    if (mtd_names[j] %in% c("G", "CG", "BG", "IBG")) {
      z_hist = mout$z
      xi_hist = mout$xi
      
      z <- get_map_labels(z_hist)$labels
      xi <- lapply(1:J, function(j) get_map_labels(sapply(xi_hist, "[[", j))$labels)
    } else {
      z <- mout$classes
      xi_j <- mout$clusters
      xi <- lapply(1:J, function(j) xi_j[[z[j]]])
    }
    
    data.frame(
      time = end_time
      , rep = rep
      , n = n
      , J = J
      , z_nmi = nett::compute_mutual_info(z, z_tru)
      , xi_nmi = hsbm::get_slice_nmi(xi, xi_tru)
      , method = mtd_names[j])
  }))
  
  out 
}, mc.cores = ncores))

res <- res %>%
  mutate(method = factor(method, labels = c("G", "CG", "BG", "IBG", "NCGE", "NCLM")))

# Visualize ----
p_z <- res %>%
  ggplot(aes(x = method, y = z_nmi, fill = method)) +
  geom_boxplot() +
  ylab(expression(bold(z)~"-NMI")) + xlab("") +
  ylim(c(0, 1)) +
  guides(fill = "none") +
  theme_minimal(base_size = 25) + theme(axis.text.x=element_text(angle = 45))

p_xi <- res %>%
  ggplot(aes(x = method, y = xi_nmi, fill = method)) +
  geom_boxplot() +
  ylab(expression(bold(xi)~"-NMI")) + xlab("") +
  ylim(c(0, 1)) +
  guides(fill = "none") +
  theme_minimal(base_size = 25) + theme(axis.text.x=element_text(angle = 45))

p_time <- res %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylab("Seconds") + xlab("") +
  guides(fill = "none") +
  scale_y_sqrt() +
  theme_minimal(base_size = 25) + theme(axis.text.x=element_text(angle = 45))

p_z + p_xi