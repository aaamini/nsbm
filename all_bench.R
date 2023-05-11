# library ---
library(igraph)
library(ggplot2)
library(dplyr)
library(parallel)
library(kableExtra)


# setwd("/project/sand/njosephs/NDP/nsbm/")
#Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
#setMethod("show", "Rcpp_NestedSBM", function(object) object$print())
source("R/nSBM_functions.R")
source("R/splice_sampler.R")
source("R/data_gen.R")
source("R/inference.R")
source("R/nsbm_wrapper.R")
source("R/NCLM.R")
source("R/setup_methods.R")

# simulation ----
#set.seed(1234)
niter = 100
K = L = 15
n_cores <- 1# detectCores()
nreps <- n_cores
rand_sbm <- TRUE
tag <- "A"

n = 100 
J = 20 
K_tru = 3 
L_tru = c(2,3,5)

set.seed(1337)

res = do.call(rbind, lapply(1:nreps, function(rep) {
  
  if (rand_sbm) {
    out = gen_rand_nsbm(n = n
                        , J = J
                        , K = K_tru
                        , L = L_tru
                        , gam = 0.4, lambda = 15)    
  } else {
    out = generate_nathans_data()  
  }
  
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
    
    if (grepl("Gibbs", mtd_names[j], fixed = TRUE)) {
      z_hist = mout$z
      xi_hist = mout$xi
      
      z <- get_map_labels(z_hist)$labels
      xi <- lapply(1:J, function(j) get_map_labels(sapply(xi_hist, "[[", j))$labels)
    } else {
      z <- mout$classes
      xi_j <- mout$clusters
      xi <- lapply(1:J, function(j) xi_j[[z[j]]])
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
}))

res$method <- factor(res$method)

if (rand_sbm) {
  state_str =  sprintf("J = %d, n = %d, nreps = %d", J, n, nreps)  
} else {
  state_str = "HSBM: Multilayer personality-friendship network"
}

res_sum <- res %>% 
  group_by(method) %>% 
  summarise(z_nmi = mean(z_nmi),
            xi_nmi = mean(xi_nmi),
            runtime = mean(runtime))

kbl(res_sum %>% arrange(desc(z_nmi)), 
                digits = 3) %>% 
  kable_paper("hover", full_width = F) %>% 
  print()