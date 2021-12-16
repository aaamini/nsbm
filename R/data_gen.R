library(nett)
library(Matrix)

symmetrize = function(A) {
  (A + t(A))/2
} 

runifmat = function(K, sym = T) {
  A = matrix(runif(K^2),K)
  if (sym) return(symmetrize(A))
}


# xi: A list.  xi[[j]] is the label vector of nodes in network j  
# z : cluster label the networks themselves.
# J : the number of networks
# L : max node community label
# K : number of network clusters
# eta: eta[[j]] connectivity matrix of network j
gen_rand_nsbm = function(n = 50, J = 10, 
                         K = 3, L = 7, 
                         lambda = 10, 
                         gam = 0.3, zeta = 0.1, 
                         pri = rep(1, L) / L, sort_z = F) {
  
  # Create connectivity matrices eta[[j]]
  eta = lapply(1:K, function(j) (1-gam)*rsymperm(L) + gam*runifmat(L))
  for (k in 1:K) {
    eta[[k]] = (1-zeta)*eta[[1]] + zeta*eta[[k]] # pull eta[[k]] towards eta[[1]]
    scale = nett::get_dcsbm_exav_deg(n, pri, eta[[k]], 1)
    eta[[k]] = pmin(eta[[k]] * lambda / scale, 1)
  }
  
  # Sample the labels, z, xi and adjacency matrices A[[j]]
  z = sample(1:K, J, replace = T)
  if (sort_z) z = sort(z)
  xi = A = vector("list", J)
  for (j in 1:J) {
    xi[[j]] = sample(L, n, replace = T, )
    A[[j]] = nett::fast_sbm(xi[[j]], eta[[z[j]]])
  }
  list(A = A, eta = eta, xi = xi, z = z)
}


generate_sparse_random_data = function(n = 50, J = 25, lambda = 10, K_tru = 3) {
  z_tru = sample(1:K_tru, J, replace = T)
  eta = lapply(1:K_tru, function(i) nett::gen_rand_conn(n, L, lambda = lambda))
  xi_tru = A = vector("list", J)
  for (j in 1:J) {
    xi_tru[[j]] = sample(L, n, replace = T, )
    A[[j]] = nett::fast_sbm(xi_tru[[j]], eta[[z_tru[j]]])
  }
  list(A = A, xi = xi_tru, z = z_tru)
}

generate_nathans_data = function(n = 40, J = 60, K = 3, sort_z = F, lambda = NULL) {
  if (K > 3) stop("K has to be <= 3.")
  eta = list(
    cbind(c(.9, .75, .5)
        , c(.75, .6, .25)
        , c(.5, .25, .1)),
    cbind(c(.1, .4, .6)
        , c(.4, .3, .1)
        , c(.6, .1, .5)),
    cbind(c(.8, .1, .3)
          , c(.1, .9, .2)
          , c(.3, .2, .7))
  )
  pri = list(c(.4, .35, .25), c(.7, .15, .15), c(.2, .4, .4))
  
  eta = eta[1:K]
  pri =pri[1:K]
  L = 3
  if (!is.null(lambda)) {
    for (k in 1:K) {
      eta[[k]] = eta[[k]] * lambda / nett::get_dcsbm_exav_deg(n, pri[[k]], eta[[k]])
    }
  }
  
  z = sample(K, J, replace = T)
  if (sort_z) z = sort(z)
  xi = A = vector("list", J)
  for (j in 1:J) {
    xi[[j]] = sample(L, n, replace = T, prob = pri[[z[j]]])
    A[[j]] = nett::fast_sbm(xi[[j]], eta[[z[j]]])
  }
  list(A = A, eta = eta, xi = xi, z = z, pri = pri)
}


# Random Graphon  ---------------------------------------------------------


# Kfun <- function(x, y, p = 2, sig2 = p /2)   exp(-sum((x-y)^2)/(2*sig2))
Kfun <- function(x, y, sig2 = 1)   exp(-mean((x-y)^2)/(2*sig2))
create_kernel_matrix = function(X, Y, sig2 = 1) {
  p = ncol(X)
  apply(Y, 1, function(y) apply(X,  1, function(x) Kfun(x, y, sig2)))
  # 1 + apply(Y, 1, function(y) apply(X,  1, function(x) Kfun(x,y)))
}


sigmoid = function(x) 1 / (1 + exp(-x))

gen_rand_graphon = function(n = 50, J = 10, 
                            K = 3, L = 7, 
                            lambda = 30, 
                            bw = 0.05, std = 2,
                            pri = lapply(1:K, function(i) rep(1, L) / L), 
                            sort_z = F) {
  
  # X = matrix(sort(runif(L)), ncol=1)
  X = matrix(seq(0,1, len = L), ncol=1)
  Kmat = create_kernel_matrix(X, X, bw)
  out = eigen(Kmat / L)  
  phi = out$vectors[,-1] # get an orthonormal basis of effectively orthogonal polynomials
  phi %>% Matrix() %>% image()
  
  eta = lapply(1:K, function(j) {
    dd = Diagonal(x = rnorm(L-1, sd = std))
    as.matrix(sigmoid( phi %*% dd %*% t(phi) ))
  })
  
  if (!is.null(lambda)) {
    for (k in 1:K) {
      eta[[k]] = eta[[k]] * lambda / nett::get_dcsbm_exav_deg(n, pri[[k]], eta[[k]])
    }
  }
  # Sample the labels, z, xi and adjacency matrices A[[j]]
  z = sample(1:K, J, replace = T)
  if (sort_z) z = sort(z)
  xi = A = vector("list", J)
  for (j in 1:J) {
    xi[[j]] = sample(L, n, replace = T, prob = pri[[z[j]]])
    A[[j]] = nett::fast_sbm(xi[[j]], eta[[z[j]]])
  }
  list(A = A, eta = eta, xi = xi, z = z)
}
 
# out = gen_rand_graphon(n = 100, J= 50, L = 5, std = 5)
# eta = out$eta
# # gr1 = create_true_fun(sig = 1)
# # u = v = seq(0,1, len = 5)
# # # outer(u, v, "*")
# # Z = matrix(gr1(expand.grid(u = u, v = v)), nrow = 5)
# library(rgl)
# # library(htmlwidgets)
# # next3d()
# persp3d(z = eta[[1]])
# eta[[1]] %>% Matrix() %>% image()
# 
# eta[[3]] %>% hist()
# norm(eta[[1]] - eta[[2]])
