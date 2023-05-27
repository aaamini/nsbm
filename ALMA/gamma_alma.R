# Libraries ----
library(parallel)
library(igraph)
library(R.matlab)

# Functions ----
source("./R/data_gen.R")

# Settings ----
ncores <- detectCores()
nreps <- 100

n <- 200          # number of nodes
J <- 20           # number of networks
K_tru <- 3        # number of true classes
L_tru <- c(2,3,5) # number of true communities in each class

# Simulation ----
runs <- expand.grid(gam = seq(0, 1, by = 0.1), rep = seq_len(nreps))

# Simulation ----
out <- mclapply(seq_len(nrow(runs)), function(ri) {
  set.seed(ri)
  
  rep <- runs[ri,"rep"]
  gam <- runs[ri, "gam"]
  
  out = gen_rand_nsbm(n = n
                      , J = J
                      , K = K_tru
                      , L = L_tru
                      , gam = gam
                      , lambda = 25)
  
  A = out$A
  z_tru = out$z
  xi_tru = out$xi
  
  list(A = A, z_tru = z_tru, xi_tru = xi_tru)
}, mc.cores = ncores)

A <- lapply(out, function(g) simplify2array(lapply(g[["A"]], as.matrix)))
A <- simplify2array(A)
A <- aperm(A, c(4, 3, 1, 2))

writeMat("./gamma_networks1.mat", A = A[1:275, , ,])
writeMat("./gamma_networks2.mat", A = A[276:550, , ,])
writeMat("./gamma_networks3.mat", A = A[551:825, , ,])
writeMat("./gamma_networks4.mat", A = A[826:1100, , ,])

z_tru <- lapply(out, "[[", 2)
xi_tru <- lapply(out, "[[", 3)

save(z_tru, xi_tru, file = "./gamma_truth.RData")
