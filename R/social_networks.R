# Libraries ----
library(igraph)
library(Matrix)
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

methods[["NCGE"]] <- methods[["NCLM"]] <- NULL
mtd_names <- names(methods)

# Data ----
K <- 5
A <- vector(mode = "list", length = K)
names(A) <- c("highschool", "facebook", "tumblr", "mit", "dblp")
for (type in names(A)) {
  
  path <- paste0("./", type, "/", type, "_ct1_")
  edges <- read.csv(paste0(path, "A.txt"), header = FALSE)
  graph_indicator <- read.csv(paste0(path, "graph_indicator.txt"), header = FALSE)
  graph_labels <- read.csv(paste0(path, "graph_labels.txt"), header = FALSE)
  
  edges <- edges %>%
    mutate(V1 = factor(V1), V2 = factor(V2))
  
  J <- length(unique(graph_indicator$V1))
  nodes <- sapply(1:J, function(j) which(graph_indicator$V1 == j))
  
  A[[type]] <- vector(mode = "list", length = J)
  for (j in 1:J) {
    ind_nodes <- which(graph_indicator$V1 == j)
    ind_edges <- which(edges$V1 %in% ind_nodes)
    G <- graph_from_edgelist(as.matrix(edges[ind_edges, ]))
    A[[type]][[j]] <- as_adj(simplify(G))
  } 
}

z_tru <- rep(1:K, times = sapply(A, length))
A <- unlist(A)

# discard the networks that break NCLM
discard <- which(is.infinite(sapply(A, function(g) sum(compute_log_moments(g)[-1]))))
A <- A[-discard]
z_tru <- z_tru[-discard]

J <- length(A)

# Cluster ----
ncores <- detectCores()
nreps <- 100
niter <- 200  # number of iteration for Gibbs samplers
K <- L <- 15  # truncation levels for NSBM models

res <- do.call(rbind, mclapply(seq_len(nreps), function(rep) {
  set.seed(rep)
  
  JJ <- sample(20:100, 1)
  ind <- sample(J, JJ)
  A <- A[ind]
  z_tru <- z_tru[ind]
  
  out = do.call(rbind, lapply(seq_along(methods), function(j) {
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    end_time = as.numeric(Sys.time() - start_time)
    
    z_hist = mout$z
    z_map <- get_map_labels(z_hist)$labels
    
    data.frame(
      time = end_time
      , rep = rep
      , J = JJ
      , iter = 1:(niter + 1)
      , z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru))
      , z_map = nett::compute_mutual_info(z_map, z_tru)
      , method = mtd_names[j])
    
  }))
  
  z_nclm <- NCLM(A)
  
  rbind(out
        , data.frame(time = NA, rep = NA, J = JJ, iter = NA
                     , z_nmi = NA
                     , z_map = nett::compute_mutual_info(z_nclm, z_tru)
                     , method = "NCLM"))
}, mc.cores = ncores))

# Visualize ----
res <- res %>%
  mutate(method = factor(method
                         , levels = c(mtd_names, "NCLM")
                         , labels = c(mtd_names, "NCLM")))

map_res <- res %>%
  group_by(method) %>%
  summarise(mean_z_map = mean(z_map)
            , lower_z = quantile(z_map, .25), upper_z = quantile(z_map, .75)
            , .groups = "drop")

mean_res =  res %>%
  filter(method != "NCLM") %>%
  group_by(iter, method) %>%
  summarise(mean_z_nmi = mean(z_nmi)
            , lower_z = quantile(z_nmi, .25), upper_z = quantile(z_nmi, .75)
            , .groups = "drop")

mean_res <- mean_res %>%
  mutate(method = factor(method
                         , levels = c(c(mtd_names, "NCLM"))
                         , labels = paste0(map_res$method, ": ", round(map_res$mean_z_map, 2), " (", round(map_res$lower_z, 2), "-", round(map_res$upper_z, 2), ")")))

mean_res %>%
  ggplot(aes(x = iter, y = mean_z_nmi, color = method)) +
  geom_line(size = 2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.25, 0.85),
    text = element_text(size = 25)
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower_z, ymax = upper_z, fill = method)
              , alpha = 0.1, linetype = "blank") +
  ylim(c(0, 1)) +  scale_x_continuous(n.breaks = 3) +
  ylab(expression(bold(z)~"-NMI")) + xlab("Iteration") +
  scale_fill_manual(values = scales::hue_pal()(7)[c(1:4, 6)], drop = FALSE) +
  scale_color_manual(values = scales::hue_pal()(7)[c(1:4, 6)], drop = FALSE)