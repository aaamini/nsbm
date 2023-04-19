# Libraries ----
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
library(patchwork)

# Functions ----
setwd("/project/sand/njosephs/NDP/nsbm/")
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

n <- 400          # number of nodes
J <- 20           # number of networks
K_tru <- 3        # number of true classes
L_tru <- c(2,3,5) # number of true communities in each class

methods[["NCGE"]] <- methods[["NCLM"]] <- NULL
mtd_names <- names(methods)

# Simulation ----
res <- do.call(rbind, mclapply(seq_len(nreps), function(rep) {
  set.seed(rep)
  
  out = gen_rand_nsbm(n = n
                      , J = J
                      , K = K_tru
                      , L = L_tru
                      , gam = 0.4
                      , lambda = 25) 
  
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
  
  out = do.call(rbind, lapply(seq_along(methods), function(j) {
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    end_time = as.numeric(Sys.time() - start_time)
    
    z_hist = mout$z
    xi_hist = mout$xi

    data.frame(
      time = end_time
      , rep = rep
      , iter = 1:(niter + 1)
      , z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru))
      , xi_nmi = sapply(xi_hist, function(xi) hsbm::get_slice_nmi(xi, xi_tru))
      , method = mtd_names[j])
  }))
  
  out 
}, mc.cores = ncores))

res <- res %>%
  mutate(method = factor(method, levels = mtd_names))

save(res, file = "./runtime_results.RData")

# Visualize ----
mean_res =  res %>%
  group_by(iter, method) %>%
  summarise(mean_z_nmi = mean(z_nmi)
            , lower_z = quantile(z_nmi, .25), upper_z = quantile(z_nmi, .75)
            , mean_xi_nmi = mean(xi_nmi)
            , lower_xi = quantile(xi_nmi, .25), upper_xi = quantile(xi_nmi, .75)
            , .groups = "drop")

p_z <- mean_res %>% 
  ggplot(aes(x = iter, y = mean_z_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.85),
    text = element_text(size = 25)
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower_z, ymax = upper_z, fill = method)
              , alpha = 0.1, linetype = "blank") +
  ylim(c(0, 1)) +
  ylab("z-NMI") + xlab("Iteration")

p_xi <- mean_res %>% 
  ggplot(aes(x = iter, y = mean_xi_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 25)) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower_xi, ymax = upper_xi, fill = method)
              , alpha = 0.1, linetype = "blank") +
  ylim(c(0, 1)) +
  ylab("xi-NMI") + xlab("Iteration")

res <- res %>%
  mutate(method = factor(method, labels = c("NCG", "CG", "BG", "IBG")))

p_time <- res %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylab("Seconds") + xlab("") +
  guides(fill = "none") +
  scale_y_sqrt() +
  theme_minimal(base_size = 25)

p_z + p_xi + p_time

ggsave("./runtime.pdf", width = 12, height = 8)
