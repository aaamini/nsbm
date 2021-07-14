comp_Bet = function(A, z, K, alpha1, beta1) {
  out = comp_blk_sums_and_sizes(A, z-1, K)  # this is from utils.cpp
  N = out$lambda
  M = out$NN
  beta(N + alpha1, M - N + beta1)
}

fit_collapsed_sbm = function(A, K, alpha1 = 1, beta1 = 1, niter = 20) {
  ut_idx = upper.tri(diag(1:K), diag = T)  # upper triangular indices
  n = nrow(A)
  z = sample(K, n, replace = T) # init labels
  z_list = vector("list", niter)
  for (it in 1:niter) {
    for (s in 1:n) {
      Bet = comp_Bet(A, z, K, alpha1, beta1)
      frq = matrix(tabulate(z, K), ncol=1)
      pri = rdirichlet(frq+1)
      
      prob = rep(0, K)
      for (rp in 1:K) {
        z_new = z
        z_new[s] = rp
        Bet_new = comp_Bet(A, z_new, K, alpha1, beta1) 
        # prob[rp] = prod((Bet_new / Bet)[ut_idx]) * (frq[rp] + 1)/(frq[z[s]])
        prob[rp] = prod((Bet_new / Bet)[ut_idx]) * pri[rp]
      } # rp
      z[s] = sample(K, 1, prob = prob) # update z
    } # s
    z_list[[it]] = list(z)
  } # it
  z_list
}
