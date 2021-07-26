add_one_to_xi = function(xi) lapply(xi, function(x) x + 1)
sub_one_from_xi = function(xi) lapply(xi, function(x) x - 1)



mat_to_list = function(x) lapply(seq_len(ncol(x)), function(i) x[,i])

splice_sampler <- function(A, K = 35, L = 55, ns, monitor = FALSE) {
  # samples from nested SBM posterior
  #
  # A         : list of adjacency matrices
  # K         : (large) number of classes
  # L         : (large) number of clusters
  # ns        : number of iterations
  # monitor   : print update every 10%
  
  require(parallel)
  n.cores <- detectCores()
  
  J <- length(A)
  n <- sapply(A, nrow)
  
  ns = ns + 1 # the first element is the initial value
  
  # pre-compute neighborhood for i = 1, ..., n_j and j = 1, ..., J
  neighbors <- lapply(A, function(A_j) apply(A_j, 1, function(u) which(u > 0)))
  
  # pre-allocate
  z <- matrix(nrow = J, ncol = ns)       # J x ns
  xi <- array(dim = c(max(n), J, ns))    # n_j x J x ns
  pi_k <- matrix(nrow = K, ncol = ns)    # K x ns
  w_k <- array(dim = c(L, K, ns))        # L x K x ns
  pi_0 <- w_0 <- log.post <- vector(length = ns)
  
  # random initialization
  z[, 1] <- sample.int(K, J, replace = TRUE)
  xi[, , 1] <- sample.int(L, J*n, replace = TRUE)
  pi_0[1] <- 1 / (J * log(J))
  w_0[1] <- max(1 / (n * log(n)))
  
  u_k <- c(rep(1/2, K - 1), 1) # mean of u_k ~ beta
  pi_k[, 1] <- sapply(seq_len(K), function(k) u_k[k] * prod(1 - u_k[seq_len(k-1)]))
  
  v_lk <- rbind(matrix(1 / (1+w_0[1]), L-1, K), 1) # mean of v_lk ~ beta
  w_k[, , 1] <- sapply(seq_len(K), function(k) {
    sapply(seq_len(L), function(l) {
      v_lk[l, k] * prod(1 - v_lk[seq_len(l-1), k])
    })})
  
  model = new(NestedSBM, A, K, L)
  model$z = z[,1] - 1
  model$xi = mat_to_list(xi[,,1]-1)
  model$w = w_k[,,1]
  model$pi = pi_k[,1]
  model$w0 = w_0[1]
  model$pi0 = pi_0[1]

  # m_hl <- lapply(seq_len(J), function(j) count_edges(A[[j]], xi[, j, 1], L))
  # log.post[1] <- ll.SBM(m_hl)
  
  # MCMC
  for (q in 2:ns) {
    
    if (monitor && q %% (ns / 10) == 0) message(q, " out of ", ns)
    
    # sample z (class assignments) ----
    z[, q] <-  z[, q-1]
    
    # p(z | pi) independent of j so repeat over J rows
    # log.pz <- matrix(log(pi_k[, q-1]), J, K, byrow = TRUE) # J x K
    
    for (j in seq_len(J)) {
      model$update_z_element_naive(j-1)
      z[j, q] = model$z[j] + 1
      # # add p(A^{j' : z_j' = k} | z, xi) and p(xi_j | w_k)
      # log.pz[j, ] <- log.pz[j, ] +
      #   sapply(seq_len(K), function(z_j) {
      #     sum(sapply(seq_len(K), function(k) {
      #       ll.SBM(m_hl[replace(z[, q], j, z_j) == k])
      #     }))}) +
      #   colSums(log(w_k[xi[seq_len(n[j]), j, q-1], , q-1]))
      
      # # sample from J different multinomials of length K
      # z[j, q] <- which(rmultinom(1, 1, exp(log.pz[j, ] - max(log.pz[j, ]))) > 0)
    }
    
    # # update number of networks in class k
    m_k <- tabulate(z[, q], nbins = K)
    
    # sample xi (cluster assignments) ----
    # log.pxi <- array(0, dim = c(L, max(n), J))
    xi[, , q] <- xi[, , q-1]

    for (k in which(m_k > 0)) {
    # cluster_update <- mclapply(which(m_k > 0), mc.cores = n.cores, function(k) { # parallel over k
      z_k <- which(z[, q] == k)
      
      # p(xi | w_lk) independent of j so repeat over L rows
      # log.pxi[, , z_k] <- log(w_k[, k, q-1]) # L x n_j
      
      for (j in z_k) {
        for (i in seq_len(n[j])) {
          model$update_xi_element(j-1, i-1)
          xi[i, j, q] = model$xi[[j]][i] +1
          # m_j <- update_edges(m_hl, neighbors, xi[seq_len(n[j]), , q], j, i)
          
          # log.pxi[, i, j] <- log.pxi[, i, j] +
          #   sapply(m_j, function(m_jl) { # p(A^{j' : z_j' = k} | z, xi)
          #     ll.SBM(replace(m_hl, j, list(m_jl))[z_k])
          #   })
          
          # # sample multinomial of length L
          # xi[i, j, q] <- which(rmultinom(1, 1
          #                                , exp(log.pxi[, i, j] - max(log.pxi[, i, j])))
          #                      > 0)
          
          # # update edge and non-edge counts for new xi
          # m_hl[[j]] <- m_j[[xi[i, j, q]]]
        } # end i loop
      } # end j loop
    #   return(list(log.pxi, xi[, , q], m_hl))
    # })
    }
       # end k loop
    
    # combine/retrieve results of mclapply
    # log.pxi <- Reduce("+", lapply(cluster_update, "[[", 1))
    
    # names(cluster_update) <- which(m_k > 0)
    # for (k in which(m_k > 0)) {
    #   clust_k <- cluster_update[[toString(k)]]
    #   for (j in which(z[, q] == k)) {
    #     xi[, j, q] <- clust_k[[2]][, j]
    #     m_hl[[j]] <- clust_k[[3]][[j]]
    #   }
    # }
    
    # sample pi_k ----
    u_k <- c(rbeta(K-1, 1 + m_k[-K]
                   , pi_0[q-1] + J - sapply(seq_len(K-1), function(k) sum(m_k[seq_len(k)])))
             , 1) # u_K = 1
    pi_k[, q] <- sapply(seq_len(K), function(k) u_k[k] * prod(1 - u_k[seq_len(k-1)]))
    
    # sample w_k ----
    n_lk <- sapply(seq_len(K), function(k) { # number of nodes in cluster l
      sapply(seq_len(L), function(l) {
        sum(xi[, z[, q] == k, q] == l)
      })})
    v_lk <- rbind(sapply(seq_len(K), function(k) {
      sapply(seq_len(L-1), function(l) {
        rbeta(1, 1 + n_lk[l, k], w_0[q-1] + sum(n_lk[(l+1):L, k]))
      })})
      , 1) # v_Lk = 1 for k = 1, ..., K
    w_k[, , q] <- sapply(seq_len(K), function(k) {
      sapply(seq_len(L), function(l) {
        v_lk[l, k] * prod(1 - v_lk[seq_len(l-1), k])
      })})
    
    # sample pi_0 ----
    pi_0[q] <- rgamma(1, K-1, -sum(log(1 - u_k[-K])))
    # pi_0[q] <- 1 / (J * log(J))
    
    # sample w_0 ----
    w_0[q] <- rgamma(1, K*(L-1), -sum(log(1 - v_lk[-L, ])))
    # w_0[q] <- max(1 / (n * log(n)))
    
    model$w = w_k[ , , q]
    model$pi = pi_k[ , q]
    model$w0 = w_0[q]
    model$pi0 = pi_0[q]
    
    # # update log post
    # log.post[q] <- sum(sapply(seq_len(J), function(j) log.pz[j, z[j, q]])) +
    #   sum(unlist(sapply(seq_len(J), function(j) diag(log.pxi[xi[seq_len(n[j]), j, q]
    #                                                          , seq_len(n[j])
    #                                                          , j]))))
    
  } # end chain loop
  
  return(list(z = z, xi = xi
              , pi_k = pi_k, w_k = w_k
              , pi_0 = pi_0, w_0 = w_0))
              # , log.post = log.post))
  
}