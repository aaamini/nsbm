# Table 3

# Libraries ----
library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)
library(parallel)
library(patchwork)
library(kableExtra)

# Functions ----
source("R/inference.R")
source("R/nsbm_wrapper.R")
source("R/NCLM.R")
source("R/setup_methods.R")

methods[["NCGE"]] <- methods[["NCLM"]] <- methods[["ALMA"]] <- NULL
mtd_names <- names(methods)

# Data ----
K <- 5
A <- vector(mode = "list", length = K)
names(A) <- c("highschool", "facebook", "tumblr", "mit", "dblp")
for (type in names(A)) {
  
  path <- paste0("./R/Applications/", type, "/", type, "_ct1_")
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
  K_tru <- length(unique(z_tru))
  
  out = do.call(rbind, lapply(seq_along(methods), function(j) {
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    end_time = as.numeric(Sys.time() - start_time)
    
    z_hist = mout$z
    z_map <- get_minVI_labels(z_hist)$labels
    
    data.frame(
      time = end_time
      , rep = rep
      , J = JJ
      , iter = 1:(niter + 1)
      , z_nmi = apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru))
      , z_map = nett::compute_mutual_info(z_map, z_tru)
      , method = mtd_names[j])
    
  }))
  
  z_nclm <- NCLM(A, K = K_tru)
  
  rbind(out
        , data.frame(time = NA, rep = NA, J = JJ, iter = NA
                     , z_nmi = NA
                     , z_map = nett::compute_mutual_info(z_nclm, z_tru)
                     , method = "NCLM"))
}, mc.cores = ncores))

# Summarize ----
map_res <- res %>%
  group_by(method) %>%
  summarise(z_nmi = mean(z_map)
            , lower_z = quantile(z_map, .25), upper_z = quantile(z_map, .75)
            , .groups = "drop")

kbl(map_res %>% arrange(desc(z_nmi)), 
    digits = 3) %>% 
  kable_paper("hover", full_width = F) %>% 
  print()
