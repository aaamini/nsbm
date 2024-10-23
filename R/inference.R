add_one_to_xi = function(xi) lapply(xi, function(x) x + 1)

get_modes <- function(x) {
  # find the modes of a vector (the values with largest frequency)
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

############### Computes the MAP estimated of the labels ##################
#' @export 
get_map_labels <- function(z, burnin=NULL){
  
  if (is.matrix(z)) {
    # single layer case --  z is an "n x niter" matrix
    n = nrow(z)
    niter = ncol(z)
    if (is.null(burnin))  burnin <- round(niter/2)

    z_cut = z[, (burnin+1):niter]
    z_map = apply(z_cut, 1, function(zrow) get_modes(zrow)[1])
    # sapply(1:n, function(i) get_modes(z_cut[i,])[1])
    z_conf = sapply(1:n, function(i) sum(z_cut[i,]==z_map[i])/length(z_cut[i,]))
    return(list(labels = z_map, confs = z_conf))
  }
  if (is.list(z)) {
    # multilayer z[[itr]][[layer_idx]]
    nlayers <- length(z[[1]])
    out = vector("list", nlayers)
    for (j in 1:nlayers) {
        z_mat <-  do.call(cbind, lapply(seq_along(z), function(itr) z[[itr]][[j]]))
        out[[j]] = get_map_labels(z_mat, burnin = burnin)
    }
    # return(out)
    return(list(labels = lapply(out, `[[`, 1), confs =  lapply(out, `[[`, 2)))
  }

  stop("z should be either a matrix or a list")
  
#   d <- list()
#   zConf <- list()
#   for (j in 1:nlayers){
#     nj <- nn[j]
#     # Zmat <- sapply((burnin+1):niter, function(itr) z[[itr]][[j]])
#     Zmat <-  do.call(cbind, lapply((burnin+1):niter, function(itr) z[[itr]][[j]])) 
#     zMAP <- sapply(1:nj, function(i)  get_modes( Zmat[i,] )[1] )
    
#     d[[j]] <- zMAP
#     zConf[[j]] <- sapply(1:nj, function(i) sum(Zmat[i,]==zMAP[i])/length(Zmat[i,]) )
    
#   }
  
#   if (consecutive) {
#     # make community labels consecutive from 1 to "number of communities"
#     dflat <- unlist(d)
#     dflat_new <- dflat
#     comm_labs <- sort(unique(dflat))
#     for (i in 1:length(comm_labs)){
#       dflat_new[dflat==comm_labs[i]] <- i
#     }
#     d <- split(dflat_new,  unlist(lapply(1:nlayers, function(i) rep(i,nn[i]))))
#   }
  
  list(labels=d, conf=zConf, Zmat=Zmat)
}

#' @export 
get_minVI_labels <- function(z, burnin=NULL){
  require(mcclust.ext)
  
  if (is.matrix(z)) {
    # single layer case --  z is an "n x niter" matrix
    n = nrow(z)
    niter = ncol(z)
    if (is.null(burnin))  burnin <- round(niter/2)
    
    z_cut = z[, (burnin+1):niter]
    z_map <- minVI(comp.psm(t(z)))$cl
    return(list(labels = z_map))
  }
  if (is.list(z)) {
    nlayers <- length(z[[1]])
    out = vector("list", nlayers)
    for (j in 1:nlayers) {
      z_mat <-  do.call(rbind, lapply(seq_along(z), function(itr) z[[itr]][[j]]))
      out[[j]] = get_minVI_labels(z_mat, burnin = burnin)
    }
    # return(out)
    return(list(labels = lapply(out, `[[`, 1)))
  }
  
  stop("z should be either a matrix or a list")
  
  list(labels=d)
}