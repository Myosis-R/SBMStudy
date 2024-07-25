# programmation de l'algorithme de latouche 2015
# utilisant l'ICL exacte pour l'inférence des attributions des groupes d'un graphe

eta_f = function(N,k,l,Z,adj_m){
  eta = 1/2+sum(matrix(rep(Z[k,],N),nrow=N)*matrix(rep(Z[l,],each=N),nrow=N)*adj_m)
  return(eta)
}

zeta_f = function(N,k,l,Z,adj_m){
  zeta = 1/2+sum(matrix(rep(Z[k,],N),nrow=N)*matrix(rep(Z[l,],each=N),nrow=N)*(1-adj_m)*(1-diag(1,N)))
  return(zeta)
}

delta_f_vec = function(k,l,g,h,i,Z,adj_m){
  delta = ((k==h)-(k==g))*Z[l,]%*%adj_m[i,]+
      ((l==h)-(l==g))*Z[k,]%*%adj_m[,i]
  return(delta)
}

rho_f_vec = function(k,l,N,g,h,i,delta,Z){
  rho = ((k==h)-(k==g))*Z[l,-i]%*%rep(1,N-1)+
    ((l==h)-(l==g))*Z[k,-i]%*%rep(1,N-1)-delta[cbind(k,l)]
  return(rho)
}


#' Implementation of the exICL algorithm for SBM
#'
#' @param N Number of nodes
#' @param Z_init A binary matrix (K,N) defining class assignment.
#' @param adj_m Adjacency matrix.
#' @param return_exicl if true returns the exact ICL value for the inference.
#' @returns A binary matrix (K,N) defining class assignment.
#' @export
exICL_f = function(N,Z_init,adj_m,return_exicl=FALSE){
  K_now = dim(Z_init)[1]
  eta = matrix(0,nrow = K_now,ncol = K_now)
  zeta = matrix(0,nrow = K_now,ncol = K_now)
  stop = 0
  Z = Z_init
  for(k in 1:K_now){
    for(l in 1:K_now){
      eta[k,l] = eta_f(N,k,l,Z,adj_m)
      zeta[k,l] = zeta_f(N,k,l,Z,adj_m)
    }
  }
  n = 1/2 + rowSums(Z)
  while (stop != 1) {
    order = sample(1:N, N, replace=F)
    stop = 1
    for(i in order){
      if(length(dim(Z))!=2) warning(Z)
      g = which.max(Z[,i])
      Delta = rep(0,K_now)
      if(n[g]<2){ # un seul sommet dans g
        for (h in seq(K_now)[-g]) {
          delta = outer(seq(K_now),seq(K_now),FUN=delta_f_vec,i=i,g=g,h=h,Z=Z,adj_m=adj_m)
          rho = outer(seq(K_now),seq(K_now),FUN=rho_f_vec,N=N,i=i,g=g,h=h,delta=delta,Z=Z)
          log_beta_1 = sum(lbeta(eta[h,-g]+delta[h,-g],zeta[h,-g]+rho[h,-g])-lbeta(eta[h,-g],zeta[h,-g]))+sum(lbeta(eta[-c(h,g),h]+delta[-c(h,g),h],zeta[-c(h,g),h]+rho[-c(h,g),h])-lbeta(eta[-c(h,g),h],zeta[-c(h,g),h]))
          log_beta_2 = 2*pi-sum(lbeta(eta[g,],zeta[g,]))-sum(lbeta(eta[-g,g],zeta[-g,g]))
          Delta[h] = log(2)+log(n[h])+lgamma((K_now-1)/2)+lgamma(K_now/2+N)-lgamma(K_now/2)-lgamma((K_now-1)/2+N)+log_beta_1+log_beta_2
        }
      }else{
        for (h in seq(K_now)[-g]) {
          delta = outer(seq(K_now),seq(K_now),FUN=delta_f_vec,i=i,g=g,h=h,Z=Z,adj_m=adj_m)
          rho = outer(seq(K_now),seq(K_now),FUN=rho_f_vec,N=N,i=i,g=g,h=h,delta=delta,Z=Z)
          log_beta_1 = sum(lbeta(eta[c(g,h),]+delta[c(g,h),],zeta[c(g,h),]+rho[c(g,h),])-lbeta(eta[c(g,h),],zeta[c(g,h),]))
          log_beta_2 = sum(lbeta(eta[-c(g,h),c(g,h)]+delta[-c(g,h),c(g,h)],zeta[-c(g,h),c(g,h)]+rho[-c(g,h),c(g,h)])-lbeta(eta[-c(g,h),c(g,h)],zeta[-c(g,h),c(g,h)]))
          Delta[h] = log(n[h]/(n[g]-1))+log_beta_1+log_beta_2
        }
      }
      idx_Delta_max = which.max(Delta)
      if(Delta[idx_Delta_max]>0){
        stop=0
        Z[g,i]=0
        Z[idx_Delta_max,i]=1
        n[g] = n[g] - 1
        n[idx_Delta_max] = n[idx_Delta_max] + 1
        if(n[g]<1){
          for(k in 1:K_now){
            eta[idx_Delta_max,k]=eta_f(N,idx_Delta_max,k,Z,adj_m)
            eta[k,idx_Delta_max]=eta_f(N,k,idx_Delta_max,Z,adj_m)
            zeta[idx_Delta_max,k]=zeta_f(N,idx_Delta_max,k,Z,adj_m)
            zeta[k,idx_Delta_max]=zeta_f(N,k,idx_Delta_max,Z,adj_m)
          }
          eta = eta[-g,-g]
          zeta = zeta[-g,-g]
          n = n[-g]
          Z = Z[-g,]
          K_now = K_now - 1
          Z = array(Z,dim=c(K_now,N))
        }else{
          for(k in 1:K_now){
            eta[idx_Delta_max,k]=eta_f(N,idx_Delta_max,k,Z,adj_m)
            eta[k,idx_Delta_max]=eta_f(N,k,idx_Delta_max,Z,adj_m)
            eta[g,k]=eta_f(N,g,k,Z,adj_m)
            eta[k,g]=eta_f(N,k,g,Z,adj_m)
            zeta[idx_Delta_max,k]=zeta_f(N,idx_Delta_max,k,Z,adj_m)
            zeta[k,idx_Delta_max]=zeta_f(N,k,idx_Delta_max,Z,adj_m)
            zeta[g,k]=zeta_f(N,g,k,Z,adj_m)
            zeta[k,g]=zeta_f(N,k,g,Z,adj_m)
          }
        }
      }
    }
  }
  if(return_exicl){
    exICL = sum(lgamma(eta)+lgamma(zeta)-lgamma(eta+zeta)-2*lgamma(1/2))+
      lgamma(K_now/2)+sum(lgamma(n))-lgamma(sum(n))-K_now*lgamma(1/2)
    return(list(Z=Z,exICL=exICL))
  }
  return(Z)
}

#' Generate random nodes' class assignment to initialize exICL algo
#'
#' @param N Number of nodes.
#' @param Ncl_init Number of nodes in each class.
#' @param adj_m Adjacency matrix.
#' @returns A binary matrix (K,N) defining class assignment.
#' @export
exICL_with_rand_init = function(N,Ncl_init,adj_m){
  Z_init = sample(rep(seq(length(Ncl_init)),times=Ncl_init))
  Z_initT = array(0,dim=c(5,N))
  Z_initT[array(c(Z_init,seq(N)),dim = c(N,2))]=1
  return(exICL_f(N,Z_initT,adj_m,TRUE))
}

#' Generate nodes' class assignment from sbm package to initialize exICL algo
#'
#' @param N Number of nodes.
#' @param adj_m Adjacency matrix.
#' @returns A binary matrix (K,N) defining class assignment.
#' @export
exICL_with_sbm_init = function(N,adj_m){
  sbm_lib <- sbm::estimateSimpleSBM(adj_m,
                                   "bernoulli",
                                   estimOptions = list(verbosity = 0, plot = FALSE))
  pi_v_init = coef(sbm_lib, 'block')
  gamma_m_init = coef(sbm_lib, 'connectivity')
  Z_init = t(sbm_lib$indMemberships)
  Z = exICL_f(N,Z_init,adj_m)
  return(Z)
}

#' For an output of SBM library calculate exICL for each model and return argmax
#' This function is badly optimize
#'
#' @param sbm output of sbm::estimateSimpleSBM
#' @returns index and value of the best model under exICL selection
#' @export
exICL_sbm = function(sbm){
  N = sbm$nbNodes
  adj_m = sbm$networkData
  nbModel = length(sbm$storedModels$indexModel)
  exICL = rep(0,nbModel)

  for(indexModel in 1:nbModel){
    K = sbm$storedModels$nbBlocks[indexModel]
    sbm$setModel(indexModel)
    Z = matrix(0,K,N)
    Z[matrix(c(sbm$memberships,1:N),nrow = N)] = 1
    eta = matrix(0,nrow = K,ncol = K)
    zeta = matrix(0,nrow = K,ncol = K)
    n = 1/2 + rowSums(Z)

    for(k in 1:K){
      for(l in 1:K){
        eta[k,l] = eta_f(N,k,l,Z,adj_m)
        zeta[k,l] = zeta_f(N,k,l,Z,adj_m)
      }
    }
    exICL[indexModel] = sum(lgamma(eta)+lgamma(zeta)-lgamma(eta+zeta)-2*lgamma(1/2))+
      lgamma(K/2)+sum(lgamma(n))-lgamma(sum(n))-K*lgamma(1/2)
  }
  idx = which.max(exICL)
  return(list(idx=idx,exICL=exICL[idx]))
}



