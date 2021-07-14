# Libraries ----
library(parallel)
library(igraph)
library(NMI)
library(mcclust.ext)

# Functions ----
source("/project/sand/njosephs/NDP/nSBM_functions.R")
source("/project/sand/njosephs/NDP/NCLM.R")

# Parameters ----
set.seed(575)
m <- seq(10, 20, 10) * 3
n <- seq(20, 100, 40)
grid <- expand.grid(m, n, stringsAsFactors = FALSE)
names(grid) <- c("m", "n")

res <- matrix(nrow = nrow(grid), ncol = 14)
colnames(res) <- c("m", "n"
                   , "nsbm_class", "ncge_class", "nclm_class"
                   , "nsbm_cluster1", "nsbm_cluster2", "nsbm_cluster3"
                   , "ncge_cluster1", "ncge_cluster2", "ncge_cluster3"
                   , "nclm_cluster1", "nclm_cluster2", "nclm_cluster3")

# (Table 1), HSBM, Paez
pm1 <- cbind(c(.9, .75, .5)
             , c(.75, .6, .25)
             , c(.5, .25, .1))
# (Section 4.2), HSBM, Paez
pm2 <- cbind(c(.8, .1, .3)
             , c(.1, .9, .2)
             , c(.3, .2, .7))
pm3 <- cbind(c(.1, .4, .6)
             , c(.4, .3, .1)
             , c(.6, .1, .5))

# Simulation ----
reps <- 1#10
for (g in 1:nrow(grid)) {
  
  m <- res[g, "m"] <- grid[g, "m"]
  n <- res[g, "n"] <- grid[g, "n"]
  
  start <- Sys.time()
  vals <- simplify2array(
    mclapply(seq_len(reps), mc.cores = detectCores(), FUN = function (r) {
      set.seed(r)
      
      G1 <- replicate(m/3, as_adj(sample_sbm(n, pm1, c(.4*n, .35*n, .25*n))))
      G2 <- replicate(m/3, as_adj(sample_sbm(n, pm2, c(.7*n, .15*n, .15*n))))
      G3 <- replicate(m/3, as_adj(sample_sbm(n, pm3, c(.2*n, .4*n, .4*n))))
      A <- c(G1, G2, G3)
      
      K.true <- c(rep(1, m/3), rep(2, m/3), rep(3, m/3))
      L.true <- rbind(
        matrix(c(rep(1, .4*n), rep(2, .35*n), rep(3, .25*n)), nrow = m/3, ncol = n, byrow = TRUE)
        , matrix(c(rep(1, .7*n), rep(2, .15*n), rep(3, .15*n)), nrow = m/3, ncol = n, byrow = TRUE)
        , matrix(c(rep(1, .2*n), rep(2, .4*n), rep(3, .4*n)), nrow = m/3, ncol = n, byrow = TRUE))
      
      # nSBM ----
      ns <- 200
      K <- L <- 10
      nsbm <- gibbs.nSBM(A, K, L, ns = ns, monitor = FALSE)
      
      # element-wise mode
      # z.mode <- apply(nsbm$z[, (.4*ns):ns], 1, function(z_j) which.max(tabulate(z_j, K)))
      # xi.mode <- sapply(seq_len(m), function(j) {
      #   sapply(seq_len(n), function (i) {
      #     which.max(tabulate(nsbm$xi[i, j, (.4*ns):ns], L))
      #   })})
      # mode via variation of information
      z.mode <- minVI(comp.psm(t(nsbm$z[, (.4*ns):ns])), t(nsbm$z[, (.4*ns):ns]), method = "draws", max.k = K)$cl
      xi.mode <- sapply(seq_len(m), function(j) {
        minVI(comp.psm(t(nsbm$xi[,j,(.4*ns):ns])), t(nsbm$xi[,j,(.4*ns):ns]), method = "draws", max.k = L)$cl
      })
      
      # NCGE -----
      ncge <- two_step(A, method = "NCGE")
      
      # NCLM -----
      nclm <- two_step(A, method = "NCLM")
      
      # NMI ----
      # class
      c(NMI(matrix(c(seq_len(m), z.mode), ncol = 2)
            , matrix(c(seq_len(m), K.true), ncol = 2))$value
        , NMI(matrix(c(seq_len(m), ncge$classes), ncol = 2)
              , matrix(c(seq_len(m), K.true), ncol = 2))$value
        , NMI(matrix(c(seq_len(m), nclm$classes), ncol = 2)
              , matrix(c(seq_len(m), K.true), ncol = 2))$value
        # communities in A1
        , mean(sapply(1:(m/3), function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), xi.mode[, j]), ncol = 2))$value
        }))
        , mean(sapply(1:(m/3), function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), ncge$clusters[[ncge$classes[j]]]), ncol = 2))$value
        }))
        , mean(sapply(1:(m/3), function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), nclm$clusters[[nclm$classes[j]]]), ncol = 2))$value
        }))
        # communities in A2
        , mean(sapply((m/3 + 1):(2*m/3), function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), xi.mode[, j]), ncol = 2))$value
        }))
        , mean(sapply((m/3 + 1):(2*m/3), function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), ncge$clusters[[ncge$classes[j]]]), ncol = 2))$value
        }))
        , mean(sapply((m/3 + 1):(2*m/3), function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), nclm$clusters[[nclm$classes[j]]]), ncol = 2))$value
        }))
        # communities in A3
        , mean(sapply((2*m/3 + 1):m, function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), xi.mode[, j]), ncol = 2))$value
        }))
        , mean(sapply((2*m/3 + 1):m, function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), ncge$clusters[[ncge$classes[j]]]), ncol = 2))$value
        }))
        , mean(sapply((2*m/3 + 1):m, function(j) {
          NMI(matrix(c(seq_len(n), L.true[j, ]), ncol = 2)
              , matrix(c(seq_len(n), nclm$clusters[[nclm$classes[j]]]), ncol = 2))$value
        })))
    }))
  
  if (is.null(dim(vals))) { # TODO: investigate errors
    message("Investigate ", g)
    vals <- simplify2array(vals[which(sapply(vals, is.numeric))])
  } else {
    message(g, " (", round(difftime(Sys.time(), start, units = "secs")), " secs) : "
            , " m = ", m
            , " , n = ", n)
  }
  
  rownames(vals) <- c("nsbm_class", "ncge_class", "nclm_class"
                      , "nsbm_cluster1", "ncge_cluster1", "nclm_cluster1"
                      , "nsbm_cluster2", "ncge_cluster2", "nclm_cluster2"
                      , "nsbm_cluster3", "ncge_cluster3", "nclm_cluster3")
  
  res[g, "nsbm_class"] <- mean(vals["nsbm_class", ])
  res[g, "ncge_class"] <- mean(vals["ncge_class", ])
  res[g, "nclm_class"] <- mean(vals["nclm_class", ])
  
  res[g, "nsbm_cluster1"] <- mean(vals["nsbm_cluster1", ])
  res[g, "nsbm_cluster2"] <- mean(vals["nsbm_cluster2", ])
  res[g, "nsbm_cluster3"] <- mean(vals["nsbm_cluster3", ])
  
  res[g, "ncge_cluster1"] <- mean(vals["ncge_cluster1", ])
  res[g, "ncge_cluster2"] <- mean(vals["ncge_cluster2", ])
  res[g, "ncge_cluster3"] <- mean(vals["ncge_cluster3", ])
  
  res[g, "nclm_cluster1"] <- mean(vals["nclm_cluster1", ])
  res[g, "nclm_cluster2"] <- mean(vals["nclm_cluster2", ])
  res[g, "nclm_cluster3"] <- mean(vals["nclm_cluster3", ])
  
}

save.image(file = "/project/sand/njosephs/NDP/nSBM_comparison.RData")


# Plot ----
library(reshape2)
library(ggplot2)

df <- melt(as.data.frame(res), id.vars = c("m", "n"))
df[c(1, 2, 4)] <- lapply(df[c(1, 2, 4)], as.numeric)
df$Method <- sub("_.*", "", df[, 3])
df$Measure <- sub(".*_", "", df[, 3])
df <- df[, c(1, 2, 4, 5, 6)]
names(df) <- c("m", "n", "NMI", "Method", "Measure")

ggplot(df, aes(x = n, y = NMI, group = Method, color = Method)) +
  ylim(c(0, 1)) +
  geom_point(aes(shape = Method)) +
  geom_line() +
  facet_grid(Measure ~ m) +
  theme_bw()
