# Libraries ----
library(igraph)
library(dplyr)
library(parallel)
library(kableExtra)

# Functions ----
source("R/inference.R")
source("R/nsbm_wrapper.R")
source("R/NCLM.R")
source("R/alma.R")
source("R/setup_methods.R")

# other methods need more than J > 2
methods[["NCLM"]] <- methods[["NCGE"]] <- methods[["ALMA"]] <- NULL
methods[["CG"]] <- NULL # some issue with (CG) when J = 2
mtd_names <- names(methods)

# Settings ----
niter <- 5000  # number of iteration for Gibbs samplers
K <- 2   # truncation levels for NSBM models (need K <= J)
L <- 15  # truncation levels for NSBM models
ncores <- detectCores()
nreps <- ncores

n <- 90           # number of nodes
J <- 2           # number of networks
K_tru <- 2        # number of true classes
L_tru <- c(3,3)   # number of true communities in each class

res = do.call(rbind, mclapply(1:nreps, function(rep) {
  set.seed(rep)
  
  eta <- list(length = 2)
  eta[[1]] <- cbind(c(.5, .1, .2)
                    , c(.1, .6, .5)
                    , c(.2, .5, .7))
  eta[[2]] <- cbind(c(.5, .1, .2)
                    , c(.1, .3, .25)
                    , c(.2, .25, .35))
  
  z_tru = rep(1:K_tru, each = ceiling(J/K_tru), length.out = J)
  xi_tru = A = vector("list", J)
  
  xi_star = lapply(seq_along(L_tru), function(k) sample(L_tru[k], n, replace=TRUE))
  for (j in 1:J) {
    xi_tru[[j]] = xi_star[[z_tru[j]]] 
    A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
  }
  
  cat(sprintf("%30s     dtime \n", "Method"))
  cat(sprintf("%30s     ----- \n", "------"))
  do.call(rbind, lapply(seq_along(methods), function(j) {
    cat(sprintf("%30s ... ", mtd_names[j]))
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    runtime = as.numeric(Sys.time() - start_time)
    
    z <- get_map_labels(mout$z)$labels
    xi <- lapply(1:J, function(j) get_map_labels(sapply(mout$xi, "[[", j))$labels)
    
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
}, mc.cores = ncores))

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
