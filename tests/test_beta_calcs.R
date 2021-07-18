Rcpp::sourceCpp("src/beta_calcs.cpp")

alpha  = 1.5
beta = 6.75
d = 3
dbar = -4 

beta(alpha + d, beta + dbar) / beta(alpha, beta)
comp_beta_ratio_v1(alpha, beta, d, dbar)
