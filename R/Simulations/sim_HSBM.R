# Figure 4

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

methods[["NCGE"]] <- methods[["ALMA"]] <- NULL
methods[["NCLM"]] <- function(A) NCLM(A)
mtd_names <- names(methods)

round_preserve_sum <- function(x, digits = 0) {
  # https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
  
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

sample_hsbm <- function(J, n, labeled = FALSE) {
  
  n <- sample(20:n, J, replace = TRUE)
  
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
    
    if (j <= J/3) {
      dist <- round_preserve_sum(c(.4*n[j], .35*n[j], .25*n[j]))
      pm <- pm1
    } else if (j <= 2*J/3) {
      dist <- round_preserve_sum(c(.7*n[j], .15*n[j], .15*n[j]))
      pm <- pm2
    } else {
      dist <- round_preserve_sum(c(.2*n[j], .4*n[j], .4*n[j]))
      pm <- pm3
    }
    
    G <- as_adj(sample_sbm(n[j], pm, dist))
    xi_tru[[j]] <- rep(1:3, times = dist)
    
    if (labeled) {
      A[[j]] <- G
    } else {
      pi <- sample(n[j])
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
n <- 100 # max number of nodes

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
      z <- get_minVI_labels(mout$z)$labels
      xi <- lapply(1:J, function(j) get_minVI_labels(sapply(mout$xi, "[[", j))$labels)
      
      xi_nmi <- hsbm::get_slice_nmi(xi, xi_tru) 
    } else {
      z <- mout
      xi_nmi <- NA
    }
    
    data.frame(
      time = end_time
      , rep = rep
      , n = n
      , J = J
      , z_nmi = nett::compute_mutual_info(z, z_tru)
      , xi_nmi = xi_nmi
      , method = mtd_names[j])
  }))
  
  out 
}, mc.cores = ncores))

# Visualize ----
res <- res %>%
  mutate(method = factor(method
                         , levels = c(sort(mtd_names[1:4]), "NCLM")
                         , labels = c(sort(mtd_names[1:4]), "NCLM")))

p_z <- res %>%
  ggplot(aes(x = method, y = z_nmi, fill = method)) +
  geom_boxplot() +
  ylab(expression(bold(z)~"-NMI")) + xlab("") +
  ylim(c(0, 1)) +
  guides(fill = "none") +
  theme_minimal(base_size = 25) + theme(axis.text.x=element_text(angle = 45)) +
  scale_fill_manual(values = scales::hue_pal()(7)[c(1:4, 7)]) +
  scale_color_manual(values = scales::hue_pal()(7)[c(1:4, 7)])

p_xi <- res %>%
  filter(method != "NCLM") %>% # nonsense for unlabeled networks
  ggplot(aes(x = method, y = xi_nmi, fill = method)) +
  geom_boxplot() +
  ylab(expression(bold(xi)~"-NMI")) + xlab("") +
  ylim(c(0, 1)) +
  guides(fill = "none") +
  theme_minimal(base_size = 25) + theme(axis.text.x=element_text(angle = 45)) +
  scale_fill_manual(values = scales::hue_pal()(7)[c(1:4, 7)]) +
  scale_color_manual(values = scales::hue_pal()(7)[c(1:4, 7)])

p_time <- res %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylab("Run Time") + xlab("") +
  guides(fill = "none") +
  theme_minimal(base_size = 25) + theme(axis.text.x=element_text(angle = 45)) +
  scale_fill_manual(values = scales::hue_pal()(7)[c(1:4, 7)]) +
  scale_color_manual(values = scales::hue_pal()(7)[c(1:4, 7)])

p_z + p_xi + p_time