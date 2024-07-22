# ensemble de fonctions utiles

#' Transform a vectorial representation of a class to a numbered representation
#'
#' @param vec A binary matrix (K,N) defining class assignment.
#' @returns A numeric vector of size N.
#' @export
#' @examples
#' \dontrun{
#' vec = matrix(c(1,0,0,1,1,0),2,3)
#' vec_to_num(vec)
#' c(1,2,1)
#' }
vec_to_num = function(vec){
  apply(vec==1,FUN=which,MARGIN=2)
}

#' Knowing nodes' class assignment and adjacency matrix, derive model parameters.
#'
#' @param Z_m A binary matrix (K,N) defining class assignment.
#' @param adj_m Adjacency matrix.
#' @returns list of SBM parameters
#' @export
Z_to_theta = function(Z_m,adj_m){
  K = nrow(Z_m)
  N = ncol(Z_m)
  gamma_m = matrix(0,K,K)
  pi_v = apply(Z_m,1,mean)
  for(k in 1:K){
    for(l in 1:K){
      temp_m = matrix(rep(Z_m[k,],N),nrow = N)*t(matrix(rep(Z_m[l,],N),nrow = N))
      gamma_m[k,l] = sum(temp_m*adj_m)/sum(temp_m*(1-diag(N)))
    }
  }
  gamma_m[is.nan(gamma_m)] = 0
  return(list(pi=pi_v,gamma=gamma_m))
}

#' Extract from colSBM output, networks' assignment.
#'
#' @param cl output of colSBM::clusterize_networks
#' @returns networks assignment
#' @export
colSBM_cluster_members <- function(cl){
  if (inherits(cl, "fitSimpleSBMPop")) return(list(cl$net_id))
  stopifnot(length(cl)==3)
  return(c(colSBM_cluster_members(cl[[2]]),
           colSBM_cluster_members(cl[[3]])))
}


colSBM_nodes_members_rec <- function(cl){
  if (inherits(cl, "fitSimpleSBMPop")) return(cl$Z)
  stopifnot(length(cl)==3)
  return(c(colSBM_nodes_members_rec(cl[[2]]),
           colSBM_nodes_members_rec(cl[[3]])))
}

#' Extract from colSBM output, nodes' assignment.
#'
#' @param cl output of colSBM::clusterize_networks
#' @returns nodes assignment
#' @export
colSBM_nodes_members <- function(cl){
  ordering_temp = unlist(colSBM_cluster_members(cl))
  ordering = rep(0,length(ordering_temp))
  for(j in 1:length(ordering_temp)){
    ordering[ordering_temp[[j]]]=j
  }
  Z_list = colSBM_nodes_members_rec(cl)
  return(Z_list[ordering])
}

#' Extract from colSBM output, networks' assignment.
#'
#' @param tree empty list
#' @param cl output of colSBM::clusterize_networks
#' @returns dendrogram
#' @export
#' @examples
#' \dontrun{
#' tree = list()
#' tree = colSBM_dendro_f(tree,cl_col)
#' class(tree) = "dendrogram"
#' plot(rev(tree))
#' }
colSBM_dendro_f <- function(tree,cl){
  if (inherits(cl, "fitSimpleSBMPop")){
    attributes(tree)=list(members=1,height=0,label=paste0(cl$net_id[1],'+'),leaf=TRUE)
  } else {
    stopifnot(length(cl)==3)
    tree[[1]] = list()
    tree[[2]] = list()
    tree[[1]] = colSBM_dendro_f(tree[[1]],cl[[2]])
    tree[[2]] = colSBM_dendro_f(tree[[2]],cl[[3]])
    cltot = cl[[1]]$BICL
    if (inherits(cl[[2]], "fitSimpleSBMPop")) cl1 = cl[[2]]$BICL else cl1 = cl[[2]][[1]]$BICL
    if (inherits(cl[[3]], "fitSimpleSBMPop")) cl2 = cl[[3]]$BICL else cl2 = cl[[3]][[1]]$BICL
    attributes(tree)=list(members=attributes(tree[[1]])$members+attributes(tree[[2]])$members,height=cltot-cl1-cl2+attributes(tree[[1]])$height+attributes(tree[[2]])$height)
  }
  return(tree)
}

#' Checking the distance between gamma and gamma with permutations
#' if the distance is low it is a sign of multiple class with the same behaviour,
#' there is a lot of limitation of this measure, especially low distance doesn't always mean harder to infer.
#'
#' @param gamma_m Probability matrix for an edge between two nodes knowing their class.
#' @returns minimun of distances
#' @export
dist_perm_gamma = function(gamma_m){
  K = nrow(gamma_m)
  if(K<2)return(0)
  perms = matrix(c(rep(1:K,each=K),rep(1:K,K)),ncol = 2)
  perms = matrix(perms[perms[,1]<perms[,2],],ncol = 2)
  linK = 1:K
  perms_m = t(apply(perms, 1, function(x,linK){linK[x]=rev(x);return(linK)},linK=linK))
  dists = apply(perms_m, 1, function(x,gamma_m){sum((gamma_m-gamma_m[x,x])**2)},gamma_m=gamma_m)
  return(min(dists))
}

#' Generate a non directed SBM
#'
#' @param N Number of nodes.
#' @param pi_v Probability vector to belong to a given class.
#' @param gamma_m Symmetric probability matrix for an edge between two nodes knowing their class.
#' @returns A binary matrix (K,N) defining class assignment and an adjacency matrix.
#' @export
rsbm_sym_f = function(N,pi_v,gamma_m){
  stopifnot(base::isSymmetric(gamma_m))
  Z = rmultinom(N,size = 1,prob = pi_v)
  Z_idx = apply(Z==1,FUN=which,MARGIN=2)
  Y = matrix(rbinom(N*N,size = 1,prob = as.vector(gamma_m[Z_idx,Z_idx])),nrow = N)
  Y = as.matrix(Matrix::tril(Y,k=-1)) + t(as.matrix(Matrix::tril(Y,k=-1)))
  return(list(Z=Z,Y=Y))
}

#' Generate a directed SBM
#'
#' @param N Number of nodes.
#' @param pi_v Probability vector to belong to a given class.
#' @param gamma_m Probability matrix for an edge between two nodes knowing their class.
#' @returns A binary matrix (K,N) defining class assignment and an adjacency matrix.
#' @export
rsbm_f = function(N,pi_v,gamma_m){
  Z = rmultinom(N,size = 1,prob = pi_v)
  Z_idx = apply(Z==1,FUN=which,MARGIN=2)
  Y = matrix(rbinom(N*N,size = 1,prob = as.vector(gamma_m[Z_idx,Z_idx])),nrow = N)
  Y[as.logical(diag(N))] = 0
  return(list(Z=Z,Y=Y))
}

