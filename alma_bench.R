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
source("R/alma_v1.R")
source("R/setup_methods2.R")

methods = methods[-c(1,2,3,5)]
mtd_names = names(methods)

# simulation ----
#set.seed(1234)
niter = 100
n_cores <- 1# detectCores()
nreps <- n_cores
rand_sbm <- TRUE
tag <- "A"

n <- 200          # number of nodes
J <- 20           # number of networks
K_tru <- 3        # number of true classes
L_tru <- c(2,3,5) # number of true communities in each class
gam <- 0.2

set.seed(1337)

res = do.call(rbind, lapply(1:nreps, function(rep) {
  
  if (rand_sbm) {
    out = gen_rand_nsbm(n = n
                        , J = J
                        , K = K_tru
                        , L = L_tru
                        , gam = gam
                        , lambda = 25)
    
  } else {
    out = generate_nathans_data()  
  }
  
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
  
  R.matlab::writeMat(paste0("./test_", tag,".mat"), 
           A = aperm(simplify2array(lapply(out[["A"]], as.matrix)), c(3, 1, 2)), 
           z_tru = z_tru, 
           xi_tru = xi_tru)
  
  cat(sprintf("%30s     dtime \n", "Method"))
  cat(sprintf("%30s     ----- \n", "------"))
  do.call(rbind, lapply(seq_along(methods), function(j) {
    cat(sprintf("%30s ... ", mtd_names[j]))
    
    start_time = Sys.time()
    mout <- methods[[j]](A)
    runtime = as.numeric(Sys.time() - start_time)
    
    cat(sprintf("%3.2f\n", runtime))
    data.frame(
      method = mtd_names[j],
      z_nmi = nett::compute_mutual_info(mout$z, z_tru),
      xi_nmi = hsbm::get_slice_nmi(mout$xi, xi_tru),
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
