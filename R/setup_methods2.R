process_two_step_labels <- function(result){
  z <- result$classes
  xi_j <- result$clusters
  xi <- lapply(1:J, function(j) xi_j[[z[j]]])
  return(list(z=z, xi=xi))
}

K = L = 10
methods = list()

methods[["G"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 3, map_labels = TRUE)
}

methods[["CG"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = T, version = 1, naive = TRUE, map_labels = TRUE)
}

methods[["BG"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 5, map_labels = TRUE)
}

methods[["IBG"]] = function(A) {
  fit_nsbm(A, K, L, niter, collapsed = F, version = 6, map_labels = TRUE)
}

methods[["NCGE"]] = function(A) {
  process_two_step_labels( two_step(A, method = "NCGE") )
}

methods[["NCLM"]] = function(A) {
  process_two_step_labels( two_step(A, method = "NCLM") )
}

methods[["ALMA"]] = function(A) {
  cat('\n')
  alma_v1(A, L_tru, verbose=TRUE)
}

methods[["ALMA (v2-paper)"]] = function(A) {
  cat('\n')
  alma(A, L_tru, init='paper', verbose=TRUE)
}

methods[["ALMA (v2-matlab)"]] = function(A) {
  cat('\n')
  alma(A, L_tru, init='matlab', verbose=TRUE)
}


mtd_names = names(methods)