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
  two_step(A, method = "NCGE")
}

methods[["NCLM"]] = function(A) {
  two_step(A, method = "NCLM")
}

mtd_names = names(methods)