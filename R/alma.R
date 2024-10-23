library(expm)
library(irlba)

# Fast projector
Kprojector <- function(X, K) {
  # projection to the nearest rank-Km matrix for each horizontal slice
  M <- dim(X)[1]
  Y <- array(0, dim = dim(X))
  
  for(i in 1:M) {
    k <- K[i]
    svdX <- irlba(X[i,,], nv = k)
    U <- svdX$u
    S <- svdX$d
    V <- svdX$v
    Y[i,,] <- U %*% diag(S) %*% t(V)
  }
  
  return(Y)
}

alma_init_matlab <- function(A, K) {
  L <- dim(A)[1]
  n <- dim(A)[2]
  M <- length(K)
  Kmax <- max(K)
  ## initialization
  
  Uhat <- apply(A, c(2, 3), sum)
  
  DEG <- apply(A, c(1, 3), sum)
  DEG <- sum(DEG)
  delta <- 2 * sqrt(Kmax) * max(DEG) / sqrt(sum(DEG^2))
  
  Uhat1 <- Uhat
  
  for(i in 1:n) {
    # d <- base::norm(Uhat[1,], type="2")
    d <- sqrt(sum(Uhat[1,]^2))
    Uhat1[i,] <- Uhat[1,] * min(delta, d) / d
  }
  
  svd1 <- svd(Uhat1)
  U <- svd1$u
  Qini <- array(0, dim = c(L, Kmax, Kmax))
  
  for(i in 1:L) {
    svd2 <- svd(t(U) %*% A[i,,] %*% U)
    U1 <- svd2$u
    Qini[i,,] <- t(U1[,1:Kmax]) %*% A[i,,] %*% U1[,1:Kmax]
  }
  
  D <- matrix(0, nrow = L, ncol = Kmax*Kmax)
  
  for(i in 1:L) {
    temp <- Qini[i,,]
    D[i,] <- c(temp)
  }
  
  return( kmeans(D, M)$cluster )
}

alma_init_paper <-  function(A, K) {
  L <- dim(A)[1]
  n <- dim(A)[2]
  
  # Assuming A is an L x n x n array
  L <- dim(A)[1]
  n <- dim(A)[2]
  M <- length(K)
  
  # Perform the 1-mode unfolding of A
  M_1 <- matrix(A, nrow = L, ncol = n^2)
  
  # Perform the truncated SVD of M_1
  svd_M_1 <- irlba(M_1, nv = M)
  
  # Store the left singular vectors 
  W0 <- svd_M_1$u
  
  return( kmeans(W0, M)$cluster )
}

orthoprojector <- function(X) {
  return(X %*% sqrtm(solve(t(X) %*% X)))
}

# K here is L_tru vector in our notation
# L here is J in our notation
alma <- function(A, K, init = 'paper', verbose = FALSE, iterMax=200, etol=1e-3) {
  if (verbose) cat('Using K = \n', K)
  
  # A is a list of L sparse n x n adjacency matrices, we turn it into a L x n x n array/tensor
  A <- simplify2array(lapply(A, as.matrix))
  A <- aperm(A, c(3, 1, 2))
  
  L <- dim(A)[1]
  n <- dim(A)[2]
  M <- length(K)
  Kmax <- max(K)
  
  if (init == "matlab") {
    if (verbose) cat('matlab init\n')
    label <- alma_init_matlab(A, K)
  } else if (init == "paper") {
    if (verbose) cat('paper init\n')
    label <- alma_init_paper(A, K)
  }
  
  Z <- nett::label_vec2mat(label)
  W <- orthoprojector(Z) # Z %*% sqrtm(solve(t(Z) %*% Z))
  
  err <- Inf
  idx <- 0
  Amatrix <- array(A, dim = c(L, n^2))
  while(err > etol && idx < iterMax) {
    idx <- idx + 1
    Wold <- W
    Qmatrix <- t(W) %*% Amatrix
    Q <- array(Qmatrix, dim = c(M, n, n))
    Q <- Kprojector(Q, K)
    Qmatrix <- array(Q, dim = c(M, n^2))
    W <- orthoprojector(Amatrix %*% t(Qmatrix))
    err <- base::norm(Wold - W, type = "F")
  }
  
  if (verbose) {
    cat('Objective error: ', err, '\n')
    cat('# of iterations: ', idx, '\n')
  }
  
  # Estimate community labels based on (Q,W)
  M <- length(K)
  z <- kmeans(W, M, iter.max = 100)$cluster
  
  label_comhat <- vector("list", M)
  for(i in seq_len(M)) {
    Qslice <- -Q[i,,]
    k <- K[i]
    svdX <- svd(Qslice)
    U <- svdX$u
    label_comhat[[i]] <- kmeans(U[,1:k], k, iter.max = 100)$cluster
  }
  return(list(z=z, xi=label_comhat[z])) # we repeat xi for all networks in the same cluster
}