library(expm)

Kprojector_v1 <- function(X, K) {
  # projectionto the nearest rank-Km matrix for each horizontal slice
  
  M <- dim(X)[1]
  Y <- array(0, dim = dim(X))
  
  for(i in 1:M) {
    k <- K[i]
    svdX <- svd(X[i,,])
    U <- svdX$u
    S <- svdX$d
    V <- svdX$v
    Y[i,,] <- U[,1:k] %*% diag(S[1:k]) %*% t(V[,1:k])
  }
  
  return(Y)
}

# K here is L_tru vector in our notation
# L here is J in our notation
AltMin_v1 <- function(A, K, verbose = FALSE) {
  # Alternating minimization
  if (verbose) cat('Using K = \n', K)
  
  # A is a list of J sparse n x n adjacency matrices, we turn it into a J x n x n array/tensor
  A <- simplify2array(lapply(A, as.matrix))
  A <- aperm(A, c(3, 1, 2))
  
  # K <- params$K
  # L <- params$L
  # M <- params$M
  # n <- params$n
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
  
  label <- kmeans(D, M)$cluster
  Z <- matrix(0, nrow = L, ncol = M)
  
  for(i in 1:L) {
    Z[i,label[i]] <- 1
  }
  
  W <- Z %*% sqrtm(solve(t(Z) %*% Z))
  
  orthoprojector <- function(X) {
    return(X %*% sqrtm(solve(t(X) %*% X)))
  }
  
  iterMax <- 200
  etol <- 1e-3
  err <- Inf
  idx <- 0
  
  Amatrix <- array(A, dim = c(L, n^2))
  
  while(err > etol && idx < iterMax) {
    idx <- idx + 1
    Wold <- W
    Qmatrix <- t(W) %*% Amatrix
    Q <- array(Qmatrix, dim = c(M, n, n))
    Q <- Kprojector_v1(Q, K)
    Qmatrix <- array(Q, dim = c(M, n^2))
    W <- orthoprojector(Amatrix %*% t(Qmatrix))
    err <- base::norm(Wold - W, type = "F")
  }
  
  if (verbose) {
    cat('Objective error: ', err, '\n')
    cat('# of iterations: ', idx, '\n')
  }
  
  return(list(Q, W))
}

alma_v1 <- function(A, K, verbose = FALSE) {
  result <- AltMin_v1(A, K, verbose = verbose)
  Q1 <- result[[1]]
  W1 <- result[[2]]
  
  M <- length(K)
  # labelhat <- kmeans(W1, M, iter.max = 100)$cluster
  z <- kmeans(W1, M, iter.max = 100)$cluster
  
  label_comhat <- vector("list", M)
  
  for(i in seq_len(M)) {
    Qslice <- -Q1[i,,]
    k <- K[i]
    svdX <- svd(Qslice)
    U <- svdX$u
    label_comhat[[i]] <- kmeans(U[,1:k], k, iter.max = 100)$cluster
  }
  return(list(z=z, xi=label_comhat[z])) # we repeat xi for all networks in the same cluster
}

