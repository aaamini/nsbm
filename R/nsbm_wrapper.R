Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

fit_nsbm = function(A, K = 10, L = 10,  
                    niter = 50, collapsed = F, naive = F, 
                    version = 1,
                    map_labels = FALSE) {
  
  model =  new(NestedSBM, A, K, L)
  
  if (collapsed) {
    if (naive) {
      fitted_model <- model$run_gibbs_naive(niter, version)
    }
    else {
      fitted_model <- model$run_gibbs(niter, version)
    }
  } else {
    fitted_model <- model$run_gibbs_via_eta(niter, version)
  }
  
  fitted_model$xi  = lapply(fitted_model$xi, function(xi) lapply(xi, function(x) x + 1))
  
  # Either output raw labels or MAP labels
  if (!map_labels) {
    return(fitted_model)
  } else{
    z_hist = fitted_model$z
    xi_hist = fitted_model$xi
    
    z <- get_map_labels(z_hist)$labels
    xi <- lapply(1:J, function(j) get_map_labels(sapply(xi_hist, "[[", j))$labels)
    return(list(z=z, xi=xi))
  }
  
}