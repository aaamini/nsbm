process_two_step_labels <- function(result){
  z <- result$classes
  xi_j <- result$clusters
  xi <- lapply(1:J, function(j) xi_j[[z[j]]])
  return(list(z=z, xi=xi))
}

methods = list()

methods[["G"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 3)
}

methods[["CG"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = T, version = 1, naive = TRUE)
}

methods[["BG"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 5)
}

methods[["IBG"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 6)
}

methods[["NCGE"]] = function(A) {
  process_two_step_labels( two_step(A, method = "NCGE", K = K_tru) )
}

methods[["NCLM"]] = function(A) {
  process_two_step_labels( two_step(A, method = "NCLM", K = K_tru) )
}

methods[["ALMA"]] = function(A) {
  cat('\n')
  alma(A, L_tru, init='matlab', verbose=FALSE)
}

mtd_names = names(methods)