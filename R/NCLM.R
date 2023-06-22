compute_log_moments <- function(A, J = 10) {
  # compute log graph moments
  # "On clustering network-valued data"
  #
  # A     : adjacency matrix
  # J     : number of moments
  
  n <- nrow(A)
  A_norm <- A / n
  
  m_k <- vector(length = J)
  A_k <- diag(n)
  
  for (k in seq_len(J)) {
    A_k <- A_k %*% A_norm
    m_k[k] <- sum(diag(A_k))  
  }
  
  return(log(m_k))
}

elbow_finder <- function(x_values, y_values) {
  # https://stackoverflow.com/a/42808962
  
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
  
  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)
  
  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }
  
  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  
  return(c(x_max_dist, y_max_dist))
}

SpecClust <- function(A, K, sphere = F, n.dim = K, lift = 0, nstart = 20) {
  
  require(mgcv)
  
  if (is.null(n.dim)) {n.dim <- K}
  n <- nrow(A)
  A <- A + matrix(lift, ncol = n, nrow = n)
  if (K == 1) {
    return(rep(1, n))
  }
  U <- slanczos(A, n.dim)$vectors
  if (sphere) {
    for (i in 1:n) {
      U[i,] = U[i,] / sqrt(sum(U[i,]^2))
    }
  }
  return(kmeans(U, K, nstart = nstart)$cluster)
}

NCGE <- function(A, K = NULL) {
  # network clustering based on graphon estimate
  # "On clustering network-valued data"
  #
  # A     : list adjacency matrix
  
  require(graphon)
  
  if (K == 1) return(rep(1, length(A)))
  
  if (!is.matrix(A[[1]])) A <- lapply(A, as.matrix)
  
  m <- length(A)
  n <- nrow(A[[1]])
  
  # 1. Graphon estimation
  P <- lapply(A, function(g) est.nbdsmooth(g)[[2]])
  
  # 2. Forming a distance matrix
  V <- sapply(P, as.vector)
  D <- as.matrix(dist(t(V))^2) / (n*(n-1)) # Frobenius distance
  
  # 3. Clustering
  if (is.null(K)) {
    kern <- exp(-D)
    lambda <- sort(abs(eigen(kern)$values), decreasing = TRUE)
    K <- elbow_finder(seq_len(m), lambda)[1] - 1
    # eigengap <- (lambda[-m] - lambda[-1]) / lambda[-1]
    # K <- which.max(eigengap)
  }
  
  return(SpecClust(D, K))
  
}

NCLM <- function(A, J = 10, K = NULL) {
  # network clustering based on log moments
  # "On clustering network-valued data"
  #
  # A     : list adjacency matrix
  # J     : max number of graph moments
  
  require(Matrix)
  
  if (K == 1) return(rep(1, length(A)))
  
  m <- length(A)
  
  # 1. Moment calculation
  m_k <- sapply(A, compute_log_moments, J)
  
  # 2. Forming a distance matrix
  D <- as.matrix(dist(t(m_k))) # Euclidean distance of log moments
  
  # 3. Clustering
  if (is.null(K)) {
    kern <- exp(-D)
    lambda <- sort(abs(eigen(kern)$values), decreasing = TRUE)
    K <- elbow_finder(seq_len(m), lambda)[1] - 1
    # eigengap <- (lambda[-m] - lambda[-1]) / lambda[-1]
    # K <- which.max(eigengap)
  }
  
  return(SpecClust(D, K))
  
}

two_step <- function(A, method, J = 10, K = NULL) {
  # cluster networks via NCLM 
  # followed by spectral clustering of average MNBS graphon within class
  #
  # A      : list adjacency matrix
  # method : NCGE or NCLM 
  # J      : max number of graph moments for NCLM
  
  require(graphon)
  
  if (!is.matrix(A[[1]])) A <- lapply(A, as.matrix)
  
  n <- nrow(A[[1]])
  
  if (method == "NCGE") {
    classes <- NCGE(A, K)
    
  } else if (method == "NCLM") {
    classes <- NCLM(A, J, K) 
  }
  
  K <- length(unique(classes))
  clusters <- vector("list", K)
  
  for (k in seq_len(K)) { # cluster within each class
    P <- lapply(A[classes == k], function(g) est.nbdsmooth(g)[[2]])
    P_avg <- Reduce("+", P) / sum(classes == k)
    
    lambda <- sort(abs(eigen(P_avg)$values), decreasing = TRUE)
    # eigengap <- (lambda[-n] - lambda[-1]) / lambda[-1]
    # ind.max <- which.max(eigengap)
    ind.max <- elbow_finder(seq_len(n), lambda)[1] - 1
    
    if (ind.max == 1) {
      clusters[[k]] <- rep(1, n)
    } else {
      clusters[[k]] <- SpecClust(A = P_avg, K = ind.max)
    }
    
  }
  
  return(list(classes = classes, clusters = clusters))
}