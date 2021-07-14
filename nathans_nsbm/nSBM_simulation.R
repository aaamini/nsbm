# library ---
library(igraph)

# functions ----
source("/project/sand/njosephs/NDP/nSBM_functions.R")

# simulation ----
set.seed(575)
m <- 60
pm1 <- cbind(c(.9, .75, .5)
             , c(.75, .6, .25)
             , c(.5, .25, .1))
pm2 <- cbind(c(.8, .1, .3)
             , c(.1, .9, .2)
             , c(.3, .2, .7))
pm3 <- cbind(c(.1, .4, .6)
             , c(.4, .3, .1)
             , c(.6, .1, .5))

K.true <- c(rep(1, m/3), rep(2, m/3), rep(3, m/3))

A <- L.true <- vector("list", length = m)

for (j in seq_len(m)) {
  
  # vary network order
  # n <- sample(seq(60, 80, 20), 1)
  n <- 80
  
  if (j <= m/3) {
    G <- as_adj(sample_sbm(n, pm1, c(.4*n, .35*n, .25*n)))
    L.true[[j]] <- c(rep(1, .4*n), rep(2, .35*n), rep(3, .25*n))
  
  } else if (j <= 2*m/3) {
    G <- as_adj(sample_sbm(n, pm2, c(.7*n, .15*n, .15*n)))
    L.true[[j]] <- c(rep(1, .7*n), rep(2, .15*n), rep(3, .15*n))
    
  } else {
    G <- as_adj(sample_sbm(n, pm3, c(.2*n, .4*n, .4*n)))
    L.true[[j]] <- c(rep(1, .2*n), rep(2, .4*n), rep(3, .4*n))
  }
  
  # unlabeled networks
  pi <- sample(n)
  A[[j]] <- G[pi, pi]
  L.true[[j]] <- L.true[[j]][pi]
}

K <- L <- 10
ns <- 250
system.time(samp <- gibbs.nSBM(A, K, L, ns = ns, monitor = TRUE))

save.image(file = "/project/sand/njosephs/NDP/nSBM_simulation.RData")

# results ----
library(NMI)
library(mcclust.ext)

# element-wise node
# z.mode <- apply(samp$z[, (.4*ns):ns], 1, function(z_j) which.max(tabulate(z_j, K)))
# table(z.mode[1:(m/3)]); table(z.mode[(m/3 + 1):(2*m/3)]); table(z.mode[(2*m/3 + 1):m])
# NMI(matrix(c(seq_len(m), z.mode), ncol = 2)
#     , matrix(c(seq_len(m), K.true), ncol = 2))$value

# mode via variation of information
z.mode <- minVI(comp.psm(t(samp$z[, (.4*ns):ns])), t(samp$z[, (.4*ns):ns]), method = "draws", max.k = K)$cl
xi.mode <- sapply(seq_len(m), function(j) {
  minVI(comp.psm(t(samp$xi[,j,(.4*ns):ns])), t(samp$xi[,j,(.4*ns):ns]), method = "draws", max.k = L)$cl
})

# plot ----
par(mfrow = c(2, 2))

nmi.z <- NMI(matrix(c(seq_len(m), K.true), ncol = 2)
             , matrix(c(seq_len(m), z.mode), ncol = 2))$value
hist(apply(samp$z[, (.4*ns):ns], 2, function(K.hat) {
    NMI(matrix(c(seq_len(m), K.true), ncol = 2)
        , matrix(c(seq_len(m), K.hat), ncol = 2))$value
  })
  , xlim = c(0, 1), xlab = "NMI"
  , main = paste0("Classes (", round(nmi.z, 3), ")"))
abline(v = nmi.z, lwd = 3)


nmi.xi1 <- sapply(seq_len(m/3), function(j) {
  NMI(matrix(c(seq_len(n), L.true[[j]]), ncol = 2)
               , matrix(c(seq_len(n), xi.mode[, j]), ncol = 2))$value
})
nmi.xi1 <- mean(nmi.xi1)
hist(simplify2array(
  mclapply(seq_len(m/3), mc.cores = detectCores(), function (j) {
    apply(samp$xi[, j, (.4*ns):ns], 2, function(L.hat) {
      NMI(matrix(c(seq_len(n), L.true[[j]]), ncol = 2)
          , matrix(c(seq_len(n), L.hat), ncol = 2))$value
    })}))
  , xlim = c(0, 1), xlab = "NMI", main = expression("Clusters in A"[1]))
abline(v = nmi.xi1, lwd = 3)

nmi.xi2 <- sapply(((m/3)+1):(2*m/3), function(j) {
  NMI(matrix(c(seq_len(n), L.true[[j]]), ncol = 2)
      , matrix(c(seq_len(n), xi.mode[, j]), ncol = 2))$value
})
nmi.xi2 <- mean(nmi.xi2)
hist(simplify2array(
  mclapply(((m/3)+1):(2*m/3), mc.cores = detectCores(), function (j) {
    apply(samp$xi[, j, (.4*ns):ns], 2, function(L.hat) {
      NMI(matrix(c(seq_len(n), L.true[[j]]), ncol = 2)
          , matrix(c(seq_len(n), L.hat), ncol = 2))$value
    })}))
  , xlim = c(0, 1), xlab = "NMI", main = expression("Clusters in A"[2]))
abline(v = nmi.xi2, lwd = 3)

nmi.xi3 <- sapply(((2*m/3)+1):m, function(j) {
  NMI(matrix(c(seq_len(n), L.true[[j]]), ncol = 2)
      , matrix(c(seq_len(n), xi.mode[, j]), ncol = 2))$value
})
nmi.xi3 <- mean(nmi.xi3)
hist(simplify2array(
  mclapply(((2*m/3)+1):m, mc.cores = detectCores(), function (j) {
    apply(samp$xi[, j, (.4*ns):ns], 2, function(L.hat) {
      NMI(matrix(c(seq_len(n), L.true[[j]]), ncol = 2)
          , matrix(c(seq_len(n), L.hat), ncol = 2))$value
    })}))
  , xlim = c(0, 1), xlab = "NMI", main = expression("Clusters in A"[3]))
abline(v = nmi.xi3, lwd = 3)
