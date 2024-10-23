# Table 3

# Libraries ----
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
library(patchwork)
# remotes::install_github("schochastics/networkdata")

# Functions ----
source("R/data_gen.R")
source("R/inference.R")
source("R/nsbm_wrapper.R")
source("R/NCLM.R")
source("R/setup_methods.R")
source("R/Applications/character_clusters.R")

methods[["NCGE"]] <- methods[["NCLM"]] <- methods[["ALMA"]] <- NULL
mtd_names <- names(methods)

# Data ----
G <- c(networkdata::starwars[1:6]
       , networkdata::got
       # , lapply(c(27, 10:11, 12:14, 28), function(i) eval(parse(text=paste0("networkdata::shakespeare_", i))))
       # , lapply(c(337:340), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Hellraiser
       # , lapply(c(367:369), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Indiana Jones
       # , lapply(c(395:397), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Jurassic Park
       , lapply(c(439:440), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # LOTR
       # , lapply(c(602:604), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Scream
       # , lapply(c(652:654), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Star Trek
       
)
J <- length(G)
z_tru <- c(rep(1, 6)  # SW
           , rep(2, 7) # GOT
           # , rep(3, 7) # Shakespeare
           , rep(4, 2) # LOTR
)
K_tru <- length(unique(z_tru))
xi_tru <- vector(mode = "list", length = 6)
for (j in 1:6) {
  g <- G[[j]]
  xi_tru[[j]] <- sw[match(V(g)$name, sw$name), "label"]
}

A <- lapply(G, as_adj)

# NSBM ----
niter <- 5000  # number of iteration for Gibbs samplers
K <- L <- 15  # truncation levels for NSBM models

set.seed(575)
res = do.call(rbind, lapply(seq_along(methods), function(j) {
  
  method <- mtd_names[j]
  
  start_time = Sys.time()
  mout <- methods[[j]](A)
  end_time = as.numeric(Sys.time() - start_time)
  
  z_hist = mout$z
  xi_hist = mout$xi
  xi_hist_known <- lapply(xi_hist, function(g) g[1:6])
  for (iter in 1:(niter + 1)) {
    for (j in 1:6) {
      ind.remove <- which(xi_tru[[j]] == 0)
      if (length(ind.remove) > 0) {
        xi_hist_known[[iter]][[j]] <- xi_hist_known[[iter]][[j]][-ind.remove]
      }
    }
  }
  for (j in 1:6) {
    ind.remove <- which(xi_tru[[j]] == 0)
    if (length(ind.remove) > 0) {
      xi_tru[[j]] <- xi_tru[[j]][-ind.remove]
    }
  }
  
  z <- get_minVI_labels(z_hist)$labels
  xi <- lapply(1:6, function(j) get_minVI_labels(sapply(xi_hist_known, "[[", j))$labels)
  
  data.frame(
    time = end_time
    , J = J
    , iter = 1:(niter + 1)
    , z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru))
    , xi_nmi = sapply(xi_hist_known, function(xi) hsbm::get_slice_nmi(xi, xi_tru))
    , z_map_nmi = nett::compute_mutual_info(z, z_tru)
    , xi_map_nmi = hsbm::get_slice_nmi(xi, xi_tru)
    , method = method)
}))

# Summarize ----
z_nclm <- NCLM(A, K = K_tru)
rbind(res %>%
        group_by(method) %>%
        summarize(mean_z = mean(z_map_nmi)
                  , mean_xi = mean(xi_map_nmi))
      , c("NCLM", nett::compute_mutual_info(z_nclm, z_tru), "NA")
)

res <- res %>%
  mutate(method = factor(method
                         , levels = c(sort(mtd_names), "NCLM")
                         , labels = c(sort(mtd_names), "NCLM")))

# Visualize ----
p_z <- res %>% 
  ggplot(aes(x = iter, y = z_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.25, 0.85),
    text = element_text(size = 25)
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylim(c(0, 1)) +  scale_x_continuous(n.breaks = 3) +
  ylab(expression(bold(z)~"-NMI")) + xlab("Iteration") +
  scale_fill_manual(values = scales::hue_pal()(7)[1:4]) +
  scale_color_manual(values = scales::hue_pal()(7)[1:4])

p_xi <- res %>%
  ggplot(aes(x = iter, y = xi_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 25)) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylim(c(0, 1)) + scale_x_continuous(n.breaks = 3) +
  ylab(expression(bold(xi)~"-NMI")) + xlab("Iteration") +scale_fill_manual(values = scales::hue_pal()(6)[1:4]) +
  scale_fill_manual(values = scales::hue_pal()(7)[1:4]) +
  scale_color_manual(values = scales::hue_pal()(7)[1:4])

p_z + p_xi