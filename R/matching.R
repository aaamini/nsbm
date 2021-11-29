get_eig_repr = function(A, K) {
  RSpectra::eigs_sym(A, K)$vectors
}
align_ortho_mats = function(X, Y) {
  svd_res = svd(t(X) %*% Y)
  svd_res$v %*% t(svd_res$u)
}

align_Y_to_X = function(X, Y) {
  Q = align_ortho_mats(X, Y)
  Y %*% Q
}

matching_perm = function(zh, z) {
  n = length(z)
  out = HungarianSolver(1-compute_confusion_matrix(zh, z)/n)
  perm = out$pairs[,2]
  perm
}

comp_all_matching_perms = function(xih, xi_tru) {
  lapply(seq_along(xih), function(j) matching_perm(xih[[j]], xi_tru[[j]]))
}

compute_entropy = function(freqs){
  probs = freqs  / sum(freqs)
  idx = probs > 0
  - sum(probs[idx]*log(probs[idx]))
}

comp_matching_score = function(xih, xi_tru, z_tru){
  ent_scores = sapply(unique(z_tru), function(k) {
    idx = z_tru == k
    list_of_perms = comp_all_matching_perms(xih[idx], xi_tru[idx])
    worst_case_ent = compute_entropy(rep(1,sum(idx)))  # entropy of the uniform
    achieved_ent = sapply(list_of_perms, paste, collapse = ", ") %>% 
      table() %>% 
      compute_entropy()  
    achieved_ent / worst_case_ent
  })
  mean(ent_scores)
  # length(unique(comp_all_matching_perms(xih, xi_tru))) # / length(xih)
}
