library(nett)

symmetrize = function(A) {
  (A + t(A))/2
} 

runifmat = function(K, sym = T) {
  A = matrix(runif(K^2),K)
  if (sym) return(symmetrize(A))
}

gen_rand_nsbm = function(n, J, K, L
                         , lambda = 10
                         , gam = 0.3
                         , labeled = TRUE) {
  # n : number of nodes per network
  # J : number of networks per class
  # K : number of classes
  # L : number of communities (in classes 1, ..., K)
  
  n.L = length(L)
  
  if (n.L == 1) {
    L <- rep(L, K)
  } else if (n.L != K) {
    stop("L needs to be a K-length vector")
  }
  
  # Create connectivity matrices eta[[j]]
  
  # eta = lapply(1:K, function(j) (1-gam)*rsymperm(L) + gam*runifmat(L))
  eta = lapply(1:K, function(j) (1-gam)*diag(L[j]) + gam*runifmat(L[j]))
  
  for (k in 1:K) {
    pri <- rep(1, L[k]) / L[k]
    scale = nett::get_dcsbm_exav_deg(n, pri, eta[[k]], 1)
    eta[[k]] = pmin(eta[[k]] * lambda / scale, 1)

  }
  
  # Sample the labels, z, xi and adjacency matrices A[[j]]
  z = rep(1:K, each = ceiling(J/K), length.out = J)
  # z = sample(1:K, J, replace = T)
  xi = A = vector("list", J)
  
  if (labeled) {
    for (j in 1:J) {
      xi[[j]] = rep(1:L[z[j]], each = ceiling(n / L[z[j]]), length.out = n)
      A[[j]] = nett::fast_sbm(xi[[j]], eta[[z[j]]])
    }
  } else {
    for (j in 1:J) {
      xi[[j]] = sample(L[z[j]], n, replace = T)
      A[[j]] = nett::fast_sbm(xi[[j]], eta[[z[j]]])
    } 
  }
  
  list(A = A, eta = eta, xi = xi, z = z)
}
