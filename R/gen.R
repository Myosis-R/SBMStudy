# fonctions de générations de SBM ou modèle de mélange de SBM
# dans l'optique de comparer des algorithmes


#' Compare different SBM inference techniques.
#'
#' @param N Number of nodes in a graph.
#' @param R Number of independent tests.
#' @returns a dataframe with multiple properties on columns and each row is a test.
#' @export
#' @importFrom foreach %dopar%
gen_sbm = function(N,R){

  K = 3

  Ncl = c(0,N%/%4,N%/%4)
  Ncl[1] = N - sum(Ncl)
  Ncl_init = rep(N%/%(K+2),K+2)
  Ncl_init[1] = N-sum(Ncl_init[-1])
  pi_v =  Ncl/N

  Ncpus <- parallel::detectCores() - 3
  par <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(par)

  iter <- itertools::isplitIndices(n=R, chunks = Ncpus)
  export <- c("dist_perm_gamma","exICL_with_rand_init","vec_to_num","exICL_f","eta_f","zeta_f","delta_f_vec","rho_f_vec")
  package <- c("sbm","graphclust")
  res_par <- foreach::foreach(I_par=iter,.export = export,.packages = package)%dopar%{

    R_par = length(I_par)

    heterogen = array(0,dim=c(R_par))
    ARI = array(0,dim=c(R_par,5))
    ICLs = array(0,dim=c(R_par,5))
    cl_diff = array(0,dim=c(R_par,5))

    for(J_par in 1:R_par){
      Z_idx =  sample(rep(c(1,2,3),times=Ncl))
      gamma_m = matrix(rbeta(K*K,1,1),nrow = K)
      adj_m = matrix(rbinom(N*N,size = 1,prob = as.vector(gamma_m[Z_idx,Z_idx])),nrow = N)
      adj_m[as.logical(diag(N))] = 0

      heterogen[J_par] = dist_perm_gamma(gamma_m)

      # sbm inference
      sbm_lib = sbm::estimateSimpleSBM(adj_m,
                                       "bernoulli",
                                       estimOptions = list(verbosity = 0, plot = FALSE))
      Z = sbm_lib$memberships
      ARI[J_par,1] = graphclust::ARI(Z,Z_idx)
      cl_diff[J_par,1] = sum(sbm_lib$blockProp>0)-sum(pi_v>0)

      # sbm inference forced with 3 classes
      sbm_lib$setModel(3)
      Z = sbm_lib$memberships
      ARI[J_par,2] = graphclust::ARI(Z,Z_idx)
      cl_diff[J_par,2] = sum(sbm_lib$blockProp>0)-sum(pi_v>0)

      # exICL inference
      exICL = exICL_with_rand_init(N,Ncl_init,adj_m)
      ARI[J_par,3] = graphclust::ARI(exICL$Z,Z_idx)
      cl_diff[J_par,3] = sum(rowSums(exICL$Z)>0)-sum(pi_v>0)
      ICLs[J_par,3] = exICL$exICL

      # exICL inference with another initialization
      exICL = exICL_with_rand_init(N,Ncl_init,adj_m)
      ARI[J_par,4] = graphclust::ARI(exICL$Z,Z_idx)
      cl_diff[J_par,4] = sum(rowSums(exICL$Z)>0)-sum(pi_v>0)
      ICLs[J_par,4] = exICL$exICL

      # exICL inference with another initialization
      exICL = exICL_with_rand_init(N,Ncl_init,adj_m)
      ARI[J_par,5] = graphclust::ARI(exICL$Z,Z_idx)
      cl_diff[J_par,5] = sum(rowSums(exICL$Z)>0)-sum(pi_v>0)
      ICLs[J_par,5] = exICL$exICL
    }
    return(list(heterogen=heterogen,ARI=ARI,cl_diff=cl_diff,ICLs=ICLs))
  }

  parallel::stopCluster(par)
  keys = names(res_par[[1]])
  dfs = lapply(keys,FUN=\(key){
    temp = purrr::map_df(res_par,\(x){dplyr::as_tibble(x[[key]],.name_repair = ~ paste0(key,1:NCOL(x[[key]])))})
    return(temp)
  })
  return(dplyr::bind_cols(dfs))
}

gen_thetaMix = function(prop,N,K){
  pi = gtools::rdirichlet(1,rep(1,K))
  while(sum(pi>0.1)<K) pi = gtools::rdirichlet(1,rep(1,K))
  gamma = matrix(rbeta(K*K,1,1),nrow = K)
  return(list(pi=pi,gamma=gamma,prop=prop))
}

rmix_sbm = function(N,labs,thetaMix){
  rsbm = lapply(labs, FUN=\(x){
    rsbm_f(N,thetaMix[[x]]$pi,thetaMix[[x]]$gamma)
  })
  return(list(listZ=lapply(rsbm, `[[`, 'Z'),listGraphs=lapply(rsbm, `[[`, 'Y'),label=labs))
}

#' Compare package graphclust and colSBM to cluster a collection of SBM.
#'
#' @param N Number of nodes in a graph.
#' @param K Maximal number of communities. (useless if thetaMix is given)
#' @param M Number of graphs.
#' @param C Number of clusters. There is, at least, one network by cluster. (useless if thetaMix is given)
#' @param R Number of independent tests.
#' @param thetaMix Optional. list of clusters properties. list(pi,gamma,prop)
#' @returns a dataframe with multiple properties on columns and each row is a test.
#' @export
#' @importFrom foreach %dopar%
gen_mix_sbm = function(N,K,M,C,R,thetaMix=NULL){

  if(! is.null(thetaMix)){
    C = length(thetaMix)
    stopifnot(M>=C)
  }

  Ncpus <- parallel::detectCores() - 3
  par <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(par)

  iter <- itertools::isplitIndices(n=R, chunks = Ncpus)
  export <- c("dist_perm_gamma","rmix_sbm","rsbm_f","gen_thetaMix","colSBM_cluster_members","colSBM_nodes_members","colSBM_nodes_members_rec")
  package <- c("sbm","graphclust","colSBM","purrr","gtools")
  res_par <- foreach::foreach(I_par=iter,.export = export,.packages = package)%dopar%{
    R_par = length(I_par)

    heterogen = array(0,dim=c(R_par))
    nodes_ARIs = array(0,dim=c(R_par,2))
    clust_ARIs = array(0,dim=c(R_par,2))
    clust_diff = array(0,dim=c(R_par,2))

    for(J_par in 1:R_par){

      if(is.null(thetaMix)){
        prop = gtools::rdirichlet(1,rep(1,C))
        thetaMix = lapply(prop,gen_thetaMix,N=N,K=K)
        heterogen[J_par] = min(sapply(thetaMix,\(x){dist_perm_gamma(x[['gamma']])}))
      } else {
        prop = lapply(thetaMix, `[[`, 'prop')
      }

      labs = c(1:C,sample(1:C,(M-C),prob=prop,replace = TRUE))

      clust = rmix_sbm(N,labs,thetaMix)

      coll_col <- colSBM::clusterize_networks(clust$listGraphs,'iid',global_opts = list(verbosity = 0,
                                                                                                       plot_details = 0,
                                                                                                       Q_max = 5))
      labs_col_temp = colSBM_cluster_members(coll_col)
      labs_col = rep(0,M)
      for(j in 1:length(labs_col_temp)){
        labs_col[labs_col_temp[[j]]]=j
      }
      clust_ARIs[J_par,1] = graphclust::ARI(clust$label,labs_col)
      clust_diff[J_par,1] = length(labs_col_temp)-C
      nodes_ARIs[J_par,1] = mean(purrr::map2_dbl(colSBM_nodes_members(coll_col),clust$listZ,graphclust::ARI))

      coll_gc <- graphclust::graphClustering(clust$listGraphs)
      clust_ARIs[J_par,2] = graphclust::ARI(clust$label,coll_gc$graphGroups)
      clust_diff[J_par,2] = length(unique(coll_gc$graphGroups))-C
      nodes_ARIs[J_par,2] = mean(map2_dbl(coll_gc$nodeClustering,clust$listZ,graphclust::ARI))

    }
    return(list(heterogen=heterogen,clust_ARIs=clust_ARIs,clust_diff=clust_diff,nodes_ARIs=nodes_ARIs))
  }

  parallel::stopCluster(par)
  keys = names(res_par[[1]])
  dfs = lapply(keys,FUN=\(key){
    temp = purrr::map_df(res_par,\(x){dplyr::as_tibble(x[[key]],.name_repair = ~ paste0(key,1:NCOL(x[[key]])))})
    return(temp)
  })
  return(dplyr::bind_cols(dfs))
}

#' Generate graphs and apply on each the vfold algorithm.
#'
#' @param N Number of nodes in a graph.
#' @param R Number of graphs generated.
#' @returns a dataframe with multiple properties on columns and each row is a test.
#' @export
#' @importFrom foreach %dopar%
gen_vfold = function(N,R){

  K = 3

  Ncpus <- parallel::detectCores() - 3
  par <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(par)

  iter <- itertools::isplitIndices(n=R, chunks = Ncpus)
  export <- c("vfold","dist_perm_gamma","rsbm_f","cross_val","method_col","method_gc","exICL_f","eta_f","zeta_f","delta_f_vec","rho_f_vec","Z_to_theta","vec_to_num")
  package <- c("sbm","gtools")
  res_par <- foreach::foreach(I_par=iter,.export = export,.packages = package)%dopar%{

    R_par = length(I_par)

    heterogen = array(0,dim=c(R_par))
    vfolds = array(0,dim=c(R_par,2))

    for(J_par in 1:R_par){
      pi_v = gtools::rdirichlet(1,rep(1,K))
      gamma_m = matrix(rbeta(K*K,1,1),nrow = K)
      rsbm = rsbm_f(N,pi_v,gamma_m)
      adj_m = rsbm[['Y']]

      heterogen[J_par] = dist_perm_gamma(gamma_m)

      vfolds[J_par,1] = vfold(adj_m,method=method_col)
      vfolds[J_par,2] = vfold(adj_m,method=method_gc)
    }
    return(list(heterogen=heterogen,vfolds=vfolds))
  }

  parallel::stopCluster(par)
  keys = names(res_par[[1]])
  dfs = lapply(keys,FUN=\(key){
    temp = purrr::map_df(res_par,\(x){dplyr::as_tibble(x[[key]],.name_repair = ~ paste0(key,1:NCOL(x[[key]])))})
    return(temp)
  })
  return(dplyr::bind_cols(dfs))
}
