# Table 3

# Libraries ----
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
library(patchwork)
# remotes::install_github("schochastics/networkdata")
library(kableExtra)

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
       , lapply(c(439:440), function(i) eval(parse(text=paste0("networkdata::movie_", i)))) # LOTR
)
J <- length(G)
z_tru <- c(rep(1, 6)  # SW
           , rep(2, 7) # GOT
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
niter <- 5000 # number of iteration for Gibbs samplers
K <- L <- 10  # truncation levels for NSBM models

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
  
  z <- get_map_labels(z_hist)$labels
  xi <- lapply(1:6, function(j) get_map_labels(sapply(xi_hist_known, "[[", j))$labels)
  
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

res_sum <- rbind(res %>%
                   group_by(method) %>%
                   summarize(z_nmi = mean(z_map_nmi)
                             , xi_nmi = mean(xi_map_nmi))
                 , c("NCLM", nett::compute_mutual_info(z_nclm, z_tru), "NA")
)
kbl(res_sum %>% arrange(desc(z_nmi)), 
    digits = 3) %>% 
  kable_paper("hover", full_width = F) %>% 
  print()
