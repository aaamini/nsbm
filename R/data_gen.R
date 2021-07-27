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
