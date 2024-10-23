# Libraries ----
library(igraph)
library(dplyr)
library(parallel)
library(kableExtra)

# Functions ----
source("R/data_gen.R")
source("R/inference.R")
source("R/nsbm_wrapper.R")
source("R/NCLM.R")
source("R/alma.R")
source("R/setup_methods.R")

# Settings ----
niter <- 100  # number of iteration for Gibbs samplers
K <- L <- 15  # truncation levels for NSBM models
ncores <- 1 #detectCores()
nreps <- ncores
labeled <- TRUE
trans_prob <- 0.7

n <- 200          # number of nodes
J <- 20           # number of networks
K_tru <- 3        # number of true classes
L_tru <- c(2,3,5) # number of true communities in each class

res = do.call(rbind, mclapply(1:nreps, function(rep) {
  set.seed(rep)
  
  out = gen_rand_nsbm(n = n
                      , J = J
                      , K = K_tru
                      , L = L_tru
                      , gam = 0.2
                      , lambda = 15
                      , labeled = labeled
                      , trans_prob = trans_prob)    
  
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
  
  cat(sprintf("%30s     dtime \n", "Method"))
  cat(sprintf("%30s     ----- \n", "------"))
  do.call(rbind, lapply(seq_along(methods), function(j) {
    cat(sprintf("%30s ... ", mtd_names[j]))
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    runtime = as.numeric(Sys.time() - start_time)
    
    if (mtd_names[j] %in% c("G", "CG", "BG", "IBG")) {
      z <- get_minVI_labels(mout$z)$labels
      xi <- lapply(1:J, function(j) get_minVI_labels(sapply(mout$xi, "[[", j))$labels)
    } else {
      z <- mout$z
      xi <- mout$xi
    }
    
    cat(sprintf("%3.2f\n", runtime))
    data.frame(
      method = mtd_names[j],
      z_nmi = nett::compute_mutual_info(z, z_tru),
      xi_nmi = hsbm::get_slice_nmi(xi, xi_tru),
      n = n,
      J = J,
      rep = rep,
      runtime = runtime
    ) 
  }))
}, mc.cores = n_cores))

res$method <- factor(res$method)

state_str =  sprintf("J = %d, n = %d, nreps = %d", J, n, nreps)  

res_sum <- res %>% 
  group_by(method) %>% 
  summarise(z_nmi = mean(z_nmi),
            xi_nmi = mean(xi_nmi),
            runtime = mean(runtime))

kbl(res_sum %>% arrange(desc(z_nmi)), 
    digits = 3) %>% 
  kable_paper("hover", full_width = F) %>% 
  print()