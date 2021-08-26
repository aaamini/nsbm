spec_net_clust = function(A, K = 3, L = 5) {

  index_seq = seq_along(A)
  xih = lapply(A, nett::spec_clust, K = L)
  Bh = lapply(index_seq, function(t) nett::estim_dcsbm(A[[t]], xih[[t]])$B)
  upt_idx = upper.tri(Bh[[1]], T)
  net_signatures = t(sapply(index_seq,  function(t) Bh[[t]][upt_idx]))
  
  list(z = kmeans(net_signatures, K)$cluster, xi = xih)
}


