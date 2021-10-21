library(nett)

symmetrize = function(A) {
  (A + t(A))/2
} 

runifmat = function(K, sym = T) {
  A = matrix(runif(K^2),K)
  if (sym) return(symmetrize(A))
}

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

generate_nathans_data = function(n = 40,J = 60) {
require(igraph)
    pm1 <- cbind(c(.9, .75, .5)
               , c(.75, .6, .25)
               , c(.5, .25, .1))
  pm2 <- cbind(c(.8, .1, .3)
               , c(.1, .9, .2)
               , c(.3, .2, .7))
  pm3 <- cbind(c(.1, .4, .6)
               , c(.4, .3, .1)
               , c(.6, .1, .5))
  
  z_tru <- c(rep(1,J/3), rep(2,J/3), rep(3,J/3))
  A <- xi_tru <- vector("list", length =J)
  for (j in seq_len(J)) {
    
    # vary network order
    # n <- sample(seq(60, 80, 20), 1)
    n <- 40
    
    if (j <=J/3) {
      G <- as_adj(sample_sbm(n, pm1, c(.4*n, .35*n, .25*n)))
     xi_tru[[j]] <- c(rep(1, .4*n), rep(2, .35*n), rep(3, .25*n))
      
    } else if (j <= 2*J/3) {
      G <- as_adj(sample_sbm(n, pm2, c(.7*n, .15*n, .15*n)))
     xi_tru[[j]] <- c(rep(1, .7*n), rep(2, .15*n), rep(3, .15*n))
      
    } else {
      G <- as_adj(sample_sbm(n, pm3, c(.2*n, .4*n, .4*n)))
     xi_tru[[j]] <- c(rep(1, .2*n), rep(2, .4*n), rep(3, .4*n))
    }
    # unlabeled networks
    pi <- sample(n)
    A[[j]] <- G[pi, pi]
   xi_tru[[j]] <-xi_tru[[j]][pi]
  }
  
  list(A = A, xi = xi_tru, z = z_tru)
}
