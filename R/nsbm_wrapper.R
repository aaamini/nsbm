Rcpp::sourceCpp("src/NestedSBM.cpp", verbose = T)
setMethod("show", "Rcpp_NestedSBM", function(object) object$print())

fit_nsbm = function(A, K = 10, L = 10,  niter = 50, collapsed = F) {
    model =  new(NestedSBM, A, K, L)
    #model$set_beta_params(1, 1)
    #model$w0 = 0.001
    #model$pi0 = 0.001
    # model
    if (collapsed) {
        # fitted_model <- model$run_gibbs_naive(niter)
        fitted_model <- model$run_gibbs(niter) 
    } else {
        fitted_model <- model$run_gibbs_via_eta(niter)
    }
    
    fitted_model$xi  = lapply(fitted_model$xi, function(xi) lapply(xi, function(x) x + 1))
    fitted_model
}

