Rcpp::sourceCpp("src/beta_calcs.cpp", verbose = T)

alpha  = 30
beta = 400
d = -3
dbar = 4 

beta(alpha + d, beta + dbar) / beta(alpha, beta)
comp_beta_ratio_v1(alpha, beta, d, dbar)
comp_beta_ratio_v2(alpha, beta, d, dbar)


alpha  = 300
beta = 4000
d = 3
dbar = -4 

beta(alpha + d, beta + dbar) / beta(alpha, beta)
comp_beta_ratio_v1(alpha, beta, d, dbar)
comp_beta_ratio_v2(alpha, beta, d, dbar)

# print_smallest_double()

alpha = 100
d = -5
comp_log_gamma_ratio(alpha, d)
log(gamma(alpha + d) / gamma(alpha))


alpha  = 300
beta = 400
d = -3
dbar = 4 

log(beta(alpha + d, beta + dbar) / beta(alpha, beta))
comp_log_beta_ratio(alpha, beta, d, dbar)


# bench::mark(v1 = beta(alpha + d, beta + dbar) / beta(alpha, beta),
#             v2 = comp_beta_ratio_v1(alpha, beta, d, dbar),
#             v3 = comp_beta_ratio_v2(alpha, beta, d, dbar))


            