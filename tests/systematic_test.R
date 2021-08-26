source("R/data_gen.R")
source("R/inference.R")
source("R/competing_methods.R")
source("R/nsbm_wrapper.R")

n = 100
K = 3
L = 5
J = 15
lambda = 15
nreps = 5
n_cores = 3
niter = 100


runs = expand.grid(rep = 1:nreps, zeta = seq(0,1,length.out=5))

methods = list()
methods[["C++ (non-collapsed)"]] = function(A) {
    out = fit_nsbm(A, K, L, niter, collapsed = F)
    list(z = get_map_labels(out$z)$labels,
        xi = get_map_labels(out$xi)$labels)
}

methods[["spec"]] = function(A) {
    spec_net_clust(A, K = 3, L = 5)
}

mtd_names = names(methods)


res = do.call(rbind, mclapply(1:nrow(runs), function(ri) {
# res = do.call(rbind, lapply(1:nrow(runs), function(ri) {
    zeta = runs[ri, "zeta"]
    rep = runs[ri, "rep"]
    out = gen_rand_nsbm(n=n, K=K, L=L, J=J, zeta=zeta, lambda=lambda)
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


state_str =  sprintf("J = %d, n = %d, nreps = %d", J, n, nreps)
p = res %>% 
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

print(p)

p = res %>% 
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

print(p)



# nett::compute_mutual_info(ztru, sp_out$z)
# # get_map_labels(out$z)
# # get_map_labels(out$xi)$labels
# # 
# mean(Matrix::rowSums(A[[1]]))
# 
# K = L = 10
# niter = 50
# out = fit_nsbm(A, niter = niter, collapsed = F)
# 
# # library(nett)
# # out$z[,niter+1]
# # model =  new(NestedSBM, A, K, L)
# # fitted_model <- model$run_gibbs(niter) 
# # nett::compute_mutual_info(ztru, fitted_model$z[,niter+1])
# zh = out$z
# plot(sapply(1:ncol(zh), function(itr) nett::compute_mutual_info(ztru, zh[,itr])))
# 
# compute_mutual_info(ztru, get_map_labels(out$z)$labels)
# 
# j = 1
# out$xi[[1]]
# 
# sourceCpp("src/dpsbm.cpp")
# aa = fit_dpsbm(A[[1]], Zcap = 2*L)
# plot(sapply(1:ncol(aa), function(itr) nett::compute_mutual_info(xi[[1]], aa[,itr])))
# compute_mutual_info(xi[[1]], get_map_labels(aa)$labels)
# 
# compute_mutual_info(xi[[1]], spec_clust(A[[1]], 7))
# 
# # samp <- splice_sampler(A, K, L, ns = niter, monitor = TRUE)
# # samp$z
