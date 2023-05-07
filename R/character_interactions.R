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

methods[["NCGE"]] <- methods[["NCLM"]] <- NULL
mtd_names <- names(methods)

# Data ----
G <- c(networkdata::starwars
       , networkdata::got
       , lapply(c(27, 10:11, 12:14, 28), function(i) eval(parse(text=paste0("networkdata::shakespeare_", i))))
       # , lapply(c(337:340), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Hellraiser
       # , lapply(c(367:369), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Indiana Jones
       # , lapply(c(395:397), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Jurassic Park
       # , lapply(c(439:440), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # LOTR
       # , lapply(c(602:604), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Scream
       # , lapply(c(652:654), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # Star Trek
       
)
J <- length(G)
z_tru <- c(rep(1, 7)  # SW
           , rep(2, 7) # GOT
           , rep(3, 7) # Shakespeare
           )

A <- lapply(G, as_adj)

# NSBM ----
niter <- 200  # number of iteration for Gibbs samplers
K <- L <- 15  # truncation levels for NSBM models

set.seed(575)
res = do.call(rbind, lapply(seq_along(methods), function(j) {
  
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
    , J = J
    , iter = 1:(niter + 1)
    , z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru))
    # , xi_nmi = sapply(xi_hist, function(xi) hsbm::get_slice_nmi(xi, xi_tru))
    # , z_map_nmi = nett::compute_mutual_info(z, z_tru)
    # , xi_map_nmi = hsbm::get_slice_nmi(xi, xi_tru)
    , method = mtd_names[j])
}))

# Summarize ----
z_nclm <- NCLM(A)
rbind(res %>%
        group_by(method) %>%
        summarize(mean_z = mean(z_nmi))
      , c("NCLM", nett::compute_mutual_info(z_nclm, z_tru))
)

res <- res %>%
  mutate(method = factor(method
                         , labels = mtd_names))

save(res, file = "./final/character_results.RData")

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
  scale_fill_manual(values = scales::hue_pal()(6)[1:4]) +
  scale_color_manual(values = scales::hue_pal()(6)[1:4])

# p_xi <- res %>% 
#   ggplot(aes(x = iter, y = xi_nmi, color = method)) +
#   geom_line(size = 1.2) +
#   theme_minimal() +
#   theme(legend.position = "none", text = element_text(size = 25)) +
#   ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
#   geom_ribbon(aes(ymin = lower_xi, ymax = upper_xi, fill = method)
#               , alpha = 0.1, linetype = "blank") +
#   ylim(c(0, 1)) + scale_x_continuous(n.breaks = 3) +
#   ylab(expression(bold(xi)~"-NMI")) + xlab("Iteration") +scale_fill_manual(values = scales::hue_pal()(6)[1:4]) +
#   scale_color_manual(values = scales::hue_pal()(6)[1:4])

p_z

ggsave("./final/character.pdf", width = 12, height = 8)
