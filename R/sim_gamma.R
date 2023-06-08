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
# source("R/setup_methods.R")
source("R/alma_v1.R")
source("R/alma_v2.R")
source("R/setup_methods3.R")


# Settings ----
niter <- 200  # number of iteration for Gibbs samplers
K <- L <- 15  # truncation levels for NSBM models
ncores <- 3 # detectCores()
nreps <- ncores # 100
labeled <- TRUE

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
                      , lambda = 25,
                      labeled = labeled)
  
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
  
  out = do.call(rbind, lapply(seq_along(methods), function(j) {
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    runtime = as.numeric(Sys.time() - start_time)
    
    # if (mtd_names[j] %in% c("G", "CG", "BG", "IBG")) {
    #   z_hist = mout$z
    #   xi_hist = mout$xi
      
    #   z <- get_map_labels(z_hist)$labels
    #   xi <- lapply(1:J, function(j) get_map_labels(sapply(xi_hist, "[[", j))$labels)
    # } else {
    #   z <- mout$classes
    #   xi_j <- mout$clusters
    #   xi <- lapply(1:J, function(j) xi_j[[z[j]]])
    # }

    data.frame(
      time = runtime,
      rep = rep,
      n = n,
      J = J,
      gam = gam,
      z_nmi = nett::compute_mutual_info(mout$z, z_tru),
      xi_nmi = hsbm::get_slice_nmi(mout$xi, xi_tru),
      method = mtd_names[j])
  }))
  
  out 
}, mc.cores = ncores))

# Load ALMA ----
# library(R.matlab)
# ALMA_res <- readMat("gamma_ALMA.mat")
# load("gamma_truth.RData")
# 
# z_ALMA <- ALMA_res[[1]]
# xi_ALMA <- ALMA_res[[2]]
# 
# ALMA_z_nmi <- sapply(1:nreps, function(rep) nett::compute_mutual_info(z_ALMA[rep, ]
#                                                                       , z_tru[[rep]]))
# ALMA_xi_nmi <- sapply(1:nreps, function(rep) hsbm::get_slice_nmi(lapply(1:J, function(j) xi_ALMA[rep,z_ALMA[rep, j],])
#                                                                  , xi_tru[[rep]]))
# 
# res <- rbind(res
#              , data.frame(time = NA
#                           , rep = 1:nreps
#                           , n = n
#                           , J = J
#                           , gam = runs$gam
#                           , z_nmi = ALMA_z_nmi
#                           , xi_nmi = ALMA_xi_nmi
#                           , method = "ALMA"))

# Visualize ----
res <- res %>%
  mutate(method = factor(method
                         , levels = mtd_names
                         , labels = mtd_names))

mean_res =  res %>% 
  group_by(method, gam) %>% 
  summarise(mean_z_nmi = mean(z_nmi)
            , lower_z = quantile(z_nmi, .25), upper_z = quantile(z_nmi, .75)
            , mean_xi_nmi = mean(xi_nmi)
            , lower_xi = quantile(xi_nmi, .25), upper_xi = quantile(xi_nmi, .75)
            , mean_time = mean(time), 
            .groups = "drop")

n_methods <- length(mtd_names)
colors <- ggsci::pal_jco()(n_methods)
p_z <- mean_res %>%
  ggplot(aes(x = gam, y = mean_z_nmi, color = method)) +
  geom_line(size = 2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.75, 0.8),
    text = element_text(size = 25)
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower_z, ymax = upper_z, fill = method)
              , alpha = 0.1, linetype = "blank") +
  ylim(c(0, 1)) +
  ylab(expression(bold(z)~"-NMI")) + xlab(expression(gamma)) 
  # scale_fill_manual(values = colors) +
  # scale_color_manual(values = colors)

p_xi <- mean_res %>% 
  ggplot(aes(x = gam, y = mean_xi_nmi, color = method)) +
  geom_line(size = 2) +
  theme_minimal() +
  theme(legend.position="none", text = element_text(size = 25)) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower_xi, ymax = upper_xi, fill = method)
              , alpha = 0.1, linetype = "blank") +
  ylim(c(0, 1)) +
  ylab(expression(bold(xi)~"-NMI")) + xlab(expression(gamma)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)

p_z + p_xi
