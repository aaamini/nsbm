spec_net_clust_naive = function(A, K = 3, L = 5) {

  index_seq = seq_along(A)
  xih = lapply(A, nett::spec_clust, K = L)
  Bh = lapply(index_seq, function(t) nett::estim_dcsbm(A[[t]], xih[[t]])$B)
  upt_idx = upper.tri(Bh[[1]], T)
  net_signatures = t(sapply(index_seq,  function(t) Bh[[t]][upt_idx]))
  
  list(z = kmeans(net_signatures, K)$cluster, xi = xih)
}

make_conn_similarity_mat = function(B_list) {
  X = do.call(rbind, lapply(B_list, as.vector))
  Di = Matrix::Matrix(as.matrix(dist(X, upper = T, diag = T)))
  Si = max(Di) - Di
  Si
}
estim_B_list = function(A, xih) {
  lapply(seq_along(A), function(t) nett::estim_dcsbm(A[[t]], xih[[t]])$B)
}

spec_net_clust = function(A, K = 3, L = 5) {
 # index_seq = seq_along(A)
  xih = lapply(A, nett::spec_clust, K = L)
  Bh_list = estim_B_list(A, xih)
  Si = make_conn_similarity_mat(Bh_list)
  
  list(z = spec_clust(Si, K), 
       xi = xih,
       Si = Si)
}



