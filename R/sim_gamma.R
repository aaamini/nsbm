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

# Settings ----
niter <- 200  # number of iteration for Gibbs samplers
K <- L <- 15  # truncation levels for NSBM models
ncores <- detectCores()
nreps <- 100

n <- 200          # number of nodes
J <- 20           # number of networks
K_tru <- 3        # number of true classes
L_tru <- c(2,3,5) # number of true communities in each class

runs <- expand.grid(gam = seq(0, 1, by = 0.1), rep = seq_len(nreps))

# Simulation ----
res <- do.call(rbind, mclapply(seq_len(nrow(runs)), function(ri) {
  set.seed(ri)
  
  rep <- runs[ri,"rep"]
  gam <- runs[ri, "gam"]
  
  out = gen_rand_nsbm(n = n
                      , J = J
                      , K = K_tru
                      , L = L_tru
                      , gam = gam
                      , lambda = 25)
  
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
      , gam = gam
      , z_nmi = nett::compute_mutual_info(z, z_tru)
      , xi_nmi = hsbm::get_slice_nmi(xi, xi_tru)
      , method = mtd_names[j])
  }))
  
  out 
}, mc.cores = ncores))

res <- res %>%
  mutate(method = factor(method, levels = mtd_names))

# Visualize ----
mean_res =  res %>% 
  group_by(method, gam) %>% 
  summarise(mean_z_nmi = mean(z_nmi)
            , lower_z = quantile(z_nmi, .25), upper_z = quantile(z_nmi, .75)
            , mean_xi_nmi = mean(xi_nmi)
            , lower_xi = quantile(xi_nmi, .25), upper_xi = quantile(xi_nmi, .75)
            , mean_time = mean(time)
            , lower_time = quantile(time, .25), upper_time = quantile(time, .75)
            , .groups = "drop")

p_z <- mean_res %>% 
  ggplot(aes(x = gam, y = mean_z_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.7, 0.85),
    text = element_text(size = 25)
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower_z, ymax = upper_z, fill = method)
              , alpha = 0.1, linetype = "blank") +
  ylim(c(0, 1)) +
  ylab(expression(bold(z)~"-NMI")) + xlab(expression(gamma))

p_xi <- mean_res %>% 
  ggplot(aes(x = gam, y = mean_xi_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  theme(legend.position="none", text = element_text(size = 25)) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower_xi, ymax = upper_xi, fill = method)
              , alpha = 0.1, linetype = "blank") +
  ylim(c(0, 1)) +
  ylab(expression(bold(xi)~"-NMI")) + xlab(expression(gamma))

p_z + p_xi