source("R/data_gen.R")
source("R/inference.R")
source("R/competing_methods.R")
source("R/nsbm_wrapper.R")
library(parallel)
library(ggplot2)
library(dplyr)

n = 100
K_tru = 2
L_tru = 5
J = 40
K = L = 10  # what we use to truncate 
lambda = 30
gam = .1
nreps = 20
n_cores = 32
niter = 100

runs = expand.grid(rep = 1:nreps, zeta = seq(0,1,length.out=8))

methods = list()
# methods[["C++ (non-collapsed)"]] = function(A) {
#     out = fit_nsbm(A, K, L, niter, collapsed = F)
#     list(z = get_map_labels(out$z)$labels,
#         xi = get_map_labels(out$xi)$labels)
# }

methods[["C++ (non-collapsed v1)"]] = function(A) {
  out = fit_nsbm(A, K, L, niter, collapsed = F, version = 1)
  list(z = get_map_labels(out$z)$labels,
       xi = get_map_labels(out$xi)$labels)
}

methods[["C++ (non-collapsed v2)"]] = function(A) {
  out = fit_nsbm(A, K, L, niter, collapsed = F, version = 2)
  list(z = get_map_labels(out$z)$labels,
       xi = get_map_labels(out$xi)$labels)
}

methods[["C++ (non-collapsed v3)"]] = function(A) {
  out = fit_nsbm(A, K, L, niter, collapsed = F, version = 3)
  list(z = get_map_labels(out$z)$labels,
       xi = get_map_labels(out$xi)$labels)
}


# Use true K and L for spec
methods[["spec"]] = function(A) {
    spec_net_clust(A, K = K_tru, L = L_tru)
}
# methods[["spec (trunc K-L)"]] = function(A) {
#   spec_net_clust(A, K = K, L = L)
# }

mtd_names = names(methods)


res = do.call(rbind, mclapply(1:nrow(runs), function(ri) {
#res = do.call(rbind, lapply(1:nrow(runs), function(ri) {
    zeta = runs[ri, "zeta"]
    rep = runs[ri, "rep"]
    out = gen_rand_nsbm(n=n, K=K_tru, L=L_tru, J=J,  lambda=lambda, gam = gam, zeta=zeta)
    A = out$A
    z_tru = out$z
    xi_tru = out$xi
    # eta_tru = out$eta
    

    do.call(rbind, lapply(seq_along(methods), function(j) { 
        dt = as.numeric(system.time( mout <- methods[[j]](A) )["elapsed"])
        data.frame(method = mtd_names[j], 
                    rep = rep,
                    zeta = zeta,
                    z_nmi = nett::compute_mutual_info(mout$z, z_tru),
                    xi_nmi = hsbm::get_slice_nmi( mout$xi, xi_tru),
                    elapsed_time = dt)
    }))
}, mc.cores = n_cores))
#}))   
ggsave(sprintf("systematic_%s_%s.png", state_str, tag), width = 10, height=5)


state_str =  sprintf("lam = %d, J = %d, n = %d, nreps = %d", 
                     lambda, J, n, nreps)

p1 = res %>% 
  group_by(zeta, method) %>% summarise(xi_nmi = mean(xi_nmi)) %>% 
  ggplot(aes(x = zeta, y = xi_nmi, color = method)) + 
  # geom_line(aes(size = method), alpha = 0.5) +
  geom_line(size = 1.2) + 
  theme_minimal() +
  # scale_colour_manual(values = c(1,1.5)) +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("xi-NMI") + xlab("zeta") + labs(title = state_str)

p2 = res %>% 
  group_by(zeta, method) %>% summarise(z_nmi = mean(z_nmi)) %>% 
  ggplot(aes(x = zeta, y = z_nmi, color = method)) + 
  # geom_line(aes(size = method), alpha = 0.5) +
  geom_line(size = 1.2) + 
  theme_minimal() +
  # scale_colour_manual(values = c(1,1.5)) +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("z-NMI") + xlab("zeta") + labs(title = state_str)


print(p1 + p2)
tag = ""
#tag = "nathan"
ggsave(sprintf("systematic_%s_%s.png", state_str, tag), width = 10, height=5)
